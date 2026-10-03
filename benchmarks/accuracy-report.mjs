// Sensitivity is visible as a moving intersection, with independent residuals.
(() => {
  const get=id=>document.getElementById(id), panel=get('accuracy-panel');
  if(!panel)return;
  const bundle=JSON.parse(get('live-data').textContent), canvas=get('accuracy-canvas'), ctx=canvas.getContext('2d');
  const controls=Array.from(panel.querySelectorAll('input,select')), fields=[];
  const presets={stable:[4,1,1,3,6,7],sensitive:[1,1,1,1.000001,2,2.000001],singular:[1,1,2,2,2,4],noisy:[1,1,1,1.000001,2,2.000101]};
  let timer, worker, workerURL, watchdog, request=0, result=null;
  const status=(text)=>{get('accuracy-status').textContent=text;};
  const format=value=>value===null?'undefined (zero baseline)':!Number.isFinite(value)?'outside float64 range':value===0?'0':Math.abs(value)<.001||Math.abs(value)>1e4?value.toExponential(3):value.toPrecision(6);
  const vector=values=>`[${values.map(format).join(', ')}]`;
  function clearResult(){
    result=null;
    for(const id of ['solution','change','condition','residual','backward','norm','rank','retained','spectrum'])get(`accuracy-${id}`).textContent='—';
    draw();
  }
  function stop(){worker?.terminate();worker=null;if(workerURL)URL.revokeObjectURL(workerURL);workerURL=null;clearTimeout(watchdog);panel.setAttribute('aria-busy','false');}
  function draw(){
    const width=Math.max(280,canvas.clientWidth||640),height=320,dpr=Math.min(2,devicePixelRatio||1);
    canvas.width=width*dpr;canvas.height=height*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,width,height);
    const style=getComputedStyle(panel), ink=style.getPropertyValue('--muted').trim()||'#75868b';
    const solutions=result?[result.base.x,...(result.perturbed?[result.perturbed.x]:[])]:[];
    const extent=Math.max(3,...solutions.flat().map(v=>Math.min(1e6,Math.abs(v)*1.25))), scale=Math.min(width,height)/2/extent*.86;
    const point=([x,y])=>[width/2+x*scale,height/2-y*scale];
    ctx.strokeStyle=ink;ctx.lineWidth=.8;ctx.beginPath();ctx.moveTo(0,height/2);ctx.lineTo(width,height/2);ctx.moveTo(width/2,0);ctx.lineTo(width/2,height);ctx.stroke();
    ctx.fillStyle=ink;ctx.font='11px ui-monospace,monospace';ctx.fillText(`x₁ / x₂ · ±${format(extent)} · auto fit`,10,18);
    if(!result)return;
    function equation(a,b,row,color,dashed){
      const u=a[row*2],v=a[row*2+1],rhs=b[row]; if(u===0&&v===0)return;
      // Clip in world coordinates, avoiding giant offscreen canvas coordinates.
      const ends=[];
      if(v!==0)for(const x of [-extent,extent]){const y=(rhs-u*x)/v;if(Number.isFinite(y)&&Math.abs(y)<=extent)ends.push([x,y]);}
      if(u!==0)for(const y of [-extent,extent]){const x=(rhs-v*y)/u;if(Number.isFinite(x)&&Math.abs(x)<=extent)ends.push([x,y]);}
      if(ends.length<2)return;
      ctx.strokeStyle=color;ctx.lineWidth=1.8;ctx.setLineDash(dashed?[5,5]:[]);ctx.beginPath();ctx.moveTo(...point(ends[0]));ctx.lineTo(...point(ends.at(-1)));ctx.stroke();ctx.setLineDash([]);
    }
    const colors=['#2a9d8f','#e58b47'];
    for(let row=0;row<2;row++){equation(result.a,result.b,row,colors[row],false);if(result.perturbed)equation(result.perturbed_a,result.perturbed_b,row,colors[row],true);}
    solutions.forEach((x,index)=>{const p=point(x);if(!p.every(Number.isFinite)||Math.abs(p[0])>width*2||Math.abs(p[1])>height*2)return;ctx.fillStyle=index?'#e25492':'#528dec';ctx.beginPath();ctx.arc(...p,index?5:7,0,Math.PI*2);ctx.fill();ctx.fillText(index?'perturbed':'baseline',p[0]+10,p[1]+(index?16:-12));});
  }
  function show(data){
    result=data.result;const r=result,b=r.base,p=r.perturbed;
    get('accuracy-condition-label').textContent=`Reciprocal condition · ${b.condition_norm||'infinity norm'}`;
    get('accuracy-norm').textContent=`${format(b.solution_norm)}${p?' → '+format(p.solution_norm):''}`;
    if(b.rank!==undefined) {
      get('accuracy-rank').textContent=`${b.rank}${p?' → '+p.rank:''}`;
      get('accuracy-retained').textContent=`${format(b.retained_reciprocal_condition)}${p?' → '+format(p.retained_reciprocal_condition):''}`;
      get('accuracy-spectrum').textContent=`${vector(b.singular_values)}${p?' → '+vector(p.singular_values):''}`;
    }
    get('accuracy-solution').textContent=p?`${vector(b.x)} → ${vector(p.x)}`:vector(b.x);
    get('accuracy-change').textContent=p?format(r.relative_solution_change):'perturbed solve rejected';
    get('accuracy-condition').textContent=`${format(b.reciprocal_condition)}${p?' → '+format(p.reciprocal_condition):''}`;
    get('accuracy-residual').textContent=`${format(b.residual_infinity)}${p?' → '+format(p.residual_infinity):''}`;
    get('accuracy-backward').textContent=`${format(b.backward_error)}${p?' → '+format(p.backward_error):''}`;
    status(`${data.checks.passed}/${data.checks.total} shared checks passed. Applied Δ = ${format(r.actual_delta)}.${r.perturbed_error?' Perturbed system: '+r.perturbed_error:''}${r.actual_delta===0&&Number(get('accuracy-delta').value)!==0?' The requested change rounded away in float64.':''}`);
    draw();
  }
  function compute(){
    clearTimeout(timer);const id=++request;clearResult();
    if(!bundle.available||!bundle.accuracy_worker_source||typeof Worker==='undefined'||typeof WebAssembly==='undefined'){
      status(bundle.reason||'Regenerate this report with a current WASM build to enable the accuracy playground.');return;
    }
    const values=fields.map(input=>input.value.trim()===''?NaN:Number(input.value));
    const input=get('accuracy-delta'),delta=input.value.trim()===''?NaN:Number(input.value);
    if(values.some(v=>!Number.isFinite(v)||Math.abs(v)>1e6)||!Number.isFinite(delta)||Math.abs(delta)>1){status('Use finite matrix entries between −10⁶ and 10⁶ and a perturbation between −1 and 1.');return;}
    try{
      if(!worker){
        workerURL=URL.createObjectURL(new Blob([bundle.accuracy_worker_source],{type:'text/javascript'}));worker=new Worker(workerURL,{type:'module',name:'matrix-accuracy'});
        worker.onmessage=({data})=>{if(data.id!==request)return;clearTimeout(watchdog);panel.setAttribute('aria-busy','false');if(data.type==='accuracy')show(data);else{clearResult();status(data.message);}};
        worker.onerror=()=>{stop();clearResult();status('The accuracy worker could not run. Reload the report to retry.');};
      }
      panel.setAttribute('aria-busy','true');status('Solving in WebAssembly…');
      clearTimeout(watchdog);watchdog=setTimeout(()=>{++request;stop();clearResult();status('Accuracy calculation timed out. Change an input to retry.');},20000);
      worker.postMessage({type:'accuracy',id,bundle,config:{a:values.slice(0,4),b:values.slice(4),delta,target:Number(get('accuracy-target').value),algorithm:get('accuracy-algorithm').value,cutoff:10**Number(get('accuracy-cutoff').value)}});
    }catch(error){stop();status(error.message);}
  }
  function cutoffControls(){
    const enabled=get('accuracy-algorithm').value==='svd';
    get('accuracy-cutoff-control').hidden=!enabled;
    get('accuracy-cutoff').disabled=!enabled;
    get('accuracy-cutoff-value').textContent=format(10**Number(get('accuracy-cutoff').value));
    panel.querySelectorAll('.accuracy-spectral').forEach(node=>node.hidden=!enabled);
  }
  function schedule(){cutoffControls();++request;clearTimeout(timer);clearTimeout(watchdog);panel.setAttribute('aria-busy','false');clearResult();status('Inputs changed; waiting for editing to pause…');timer=setTimeout(compute,300);}
  for(let i=0;i<6;i++){
    const input=document.createElement('input');input.type='number';input.step='any';input.min='-1000000';input.max='1000000';input.setAttribute('aria-label',i<4?`A row ${Math.floor(i/2)+1} column ${i%2+1}`:`b row ${i-3}`);input.addEventListener('input',()=>{get('accuracy-preset').value='custom';schedule();});fields.push(input);get(i<4?'accuracy-a':'accuracy-b').append(input);
  }
  function preset(){const values=presets[get('accuracy-preset').value];if(!values)return;fields.forEach((input,i)=>input.value=values[i]);if(get('accuracy-preset').value==='noisy')get('accuracy-algorithm').value='svd';schedule();}
  get('accuracy-preset').addEventListener('change',preset);
  controls.filter(input=>input.id!=='accuracy-preset').forEach(input=>input.addEventListener(input.tagName==='SELECT'?'change':'input',schedule));
  addEventListener('resize',draw);addEventListener('pagehide',()=>{clearTimeout(timer);stop();});
  new MutationObserver(draw).observe(document.documentElement,{attributes:true,attributeFilter:['data-theme']});
  preset();
})();
