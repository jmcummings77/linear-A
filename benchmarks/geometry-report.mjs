// Interactive geometry from verified WASM outputs and optional captured C steps.
(() => {
  const get = id => document.getElementById(id);
  const bundle = JSON.parse(get('live-data').textContent);
  const canvas = get('geometry-canvas'), ctx = canvas.getContext('2d');
  const identity = [1,0,0, 0,1,0, 0,0,1];
  const rotation = [0,-1,0, 1,0,0, 0,0,1], shear = [1,1,0, 0,1,0, 0,0,1];
  const presets = {
    'rotate-shear': [rotation, shear], 'shear-rotate': [shear, rotation],
    reflection: [identity, [-1,0,0, 0,1,0, 0,0,1]],
    collapse: [identity, [1,1,0, 0,0,0, 0,0,1]],
    'symmetric-stretch': [[2,1,0, 1,2,0, 0,0,.5], [1,0,0, 0,1,0, 0,0,3]],
    'axis-rotation': [identity, identity], 'cross-direction': [identity, identity],
  };
  const symbols = {multiply:'MN',add:'M + N',subtract:'M − N',transpose:'Mᵀ',scale:'sM'};
  const fields = {m:[],n:[]}, cells = [], sourceLines = new Map();
  const visibility = {axes:get('geometry-show-axes'),box:get('geometry-show-box'),ellipse:get('geometry-show-ellipse'),grid:get('geometry-show-grid')};
  const appearance = {style:get('geometry-arrow-style'),color:get('geometry-arrow-color'),motion:get('geometry-arrow-motion'),width:get('geometry-arrow-width')};
  const autoOrbit = get('geometry-auto-orbit');
  const reducedMotion = matchMedia('(prefers-reduced-motion: reduce)');
  const arrowColors = {cyan:[104,225,232],blue:[138,184,248],coral:[255,151,139],violet:[209,164,255],white:[235,243,244]};
  const axisInputs = Array.from(get('geometry-axis').querySelectorAll('input'));
  let scene = null, worker = null, workerUrl = null, timeout = null;
  let progress = 0, step = 0, playing = false, animation = null, previousTime = null;
  let arrowPhase = 0, arrowPreviousTime = null;
  let orbitPreviousTime = null;
  let yaw = -.68, pitch = .43, zoom = 1, drag = null;
  const minZoom = .45, maxZoom = 25;
  const changeZoom = factor => {zoom=Math.max(minZoom,Math.min(maxZoom,zoom*factor));render();};
  const fmt = value => Math.abs(value) >= 1e4 || (value !== 0 && Math.abs(value) < 1e-4)
    ? value.toExponential(2) : Number(value.toFixed(4)).toString();
  const determinant = a => a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6]);
  const transform = (a,v) => [0,1,2].map(r=>a[r*3]*v[0]+a[r*3+1]*v[1]+a[r*3+2]*v[2]);
  const add = (a,b) => a.map((x,i)=>x+b[i]);
  const scale = (v,s) => v.map(x=>x*s);
  const isSteps = () => get('geometry-mode').value === 'steps';
  const status = (message,state='') => {get('geometry-status').textContent=message;get('geometry-status').className=`small geometry-status ${state}`;};

  for (const name of ['m','n']) for (let i=0;i<9;i++) {
    const input=document.createElement('input');
    input.type='number'; input.min='-10'; input.max='10'; input.step='0.25';
    input.setAttribute('aria-label',`${name.toUpperCase()} row ${Math.floor(i/3)+1} column ${i%3+1}`);
    input.addEventListener('input',()=>{get('geometry-preset').value='custom';invalidate();});
    fields[name].push(input);get(`geometry-${name}`).append(input);
  }
  for(let i=0;i<9;i++){const cell=document.createElement('span');cell.textContent='—';get('geometry-matrix').append(cell);cells.push(cell);}
  if(bundle.geometry_source){
    get('geometry-source-title').textContent=bundle.geometry_source.path+' · current trace build';
    bundle.geometry_source.lines.forEach((line,i)=>{
      const number=bundle.geometry_source.start_line+i, row=document.createElement('span'),label=document.createElement('b');
      label.textContent=number;row.append(label,document.createTextNode(line));get('geometry-source').append(row);sourceLines.set(number,row);
    });
  }

  const arrowsMoving = () => appearance.motion.value !== 'static' && !reducedMotion.matches;
  const orbitMoving = () => autoOrbit.checked && !reducedMotion.matches && !drag;
  function syncAnimation(){
    const active = scene && !document.hidden && (playing || arrowsMoving() || orbitMoving());
    if(active && animation === null) animation=requestAnimationFrame(tick);
    if(!active){cancelAnimationFrame(animation);animation=null;arrowPreviousTime=null;orbitPreviousTime=null;}
  }
  function stopAnimation(){playing=false;previousTime=null;syncAnimation();get('geometry-play').textContent=isSteps()?'Play calculation':'Play transformation';}
  function endWorker(){worker?.terminate();worker=null;if(workerUrl)URL.revokeObjectURL(workerUrl);workerUrl=null;clearTimeout(timeout);get('geometry-compute').textContent='Compute in WebAssembly';}
  function invalidate(){endWorker();stopAnimation();scene=null;progress=0;step=0;get('geometry-play').disabled=true;get('geometry-progress').disabled=true;get('geometry-mode').value='morph';get('geometry-mode').options[1].disabled=true;get('geometry-step-controls').hidden=true;get('geometry-progress').parentElement.hidden=false;get('geometry-badge').textContent='INPUTS CHANGED';status('Compute to apply these inputs.');render();}
  function applyPreset(){
    const preset=get('geometry-preset').value,pair=presets[preset];if(!pair)return;
    pair.forEach((values,j)=>values.forEach((v,i)=>{fields[j?'n':'m'][i].value=v;}));
    get('geometry-operation').value=preset==='axis-rotation'?'rotate':preset==='cross-direction'?'cross':'multiply';
    if(preset==='axis-rotation'||preset==='cross-direction'){
      axisInputs.forEach((input,i)=>{input.value=i===2?'1':'0';});get('geometry-angle').value='90';
    }
    operationChanged();
  }
  function operationChanged(){
    const op=get('geometry-operation').value;
    fields.n.forEach(input=>{input.disabled=!['multiply','add','subtract'].includes(op);});
    get('geometry-scalar-label').hidden=op!=='scale';
    get('geometry-vector-controls').hidden=!['rotate','cross'].includes(op);
    get('geometry-angle-label').hidden=op!=='rotate';
    get('geometry-vector-note').textContent=op==='rotate'?'Right-hand rule about u. Its length is ignored; its direction must be nonzero.':
      'Each field vector Mx is crossed with u: (Mx) × u. The length of u scales the result. This is not a cross product of two arbitrary matrices.';
    invalidate();
  }
  function configuration(){
    const read=inputs=>inputs.map(input=>{if(input.value.trim()==='')throw new Error('Every matrix cell needs a number.');const n=Number(input.value);if(!Number.isFinite(n)||Math.abs(n)>10)throw new Error('Use finite matrix entries between −10 and 10.');return n;});
    const operation=get('geometry-operation').value,scalar=operation==='scale'?Number(get('geometry-scalar').value):1.25;
    if((operation==='scale'&&get('geometry-scalar').value.trim()==='')||!Number.isFinite(scalar)||Math.abs(scalar)>10)throw new Error('Use a scalar between −10 and 10.');
    const config={operation,m:read(fields.m),n:['multiply','add','subtract'].includes(operation)?read(fields.n):identity.slice(),scalar};
    if(operation==='rotate'||operation==='cross')config.axis=read(axisInputs);
    if(operation==='rotate'){
      const input=get('geometry-angle'),degrees=Number(input.value);
      if(input.value.trim()===''||!Number.isFinite(degrees)||Math.abs(degrees)>360)throw new Error('Use an angle from −360 to 360 degrees.');
      if(!config.axis.some(x=>x!==0))throw new Error('Rotation axis must be nonzero.');
      config.radians=degrees*Math.PI/180;
    }
    return config;
  }
  function compute(){
    if(worker){endWorker();status('Computation cancelled.');return;}
    let config;try{config=configuration();}catch(error){status(error.message,'failed');return;}
    stopAnimation();scene=null;render();status('Loading WASM and verifying the shared math cases…');
    try{
      workerUrl=URL.createObjectURL(new Blob([bundle.geometry_worker_source],{type:'text/javascript'}));
      const current=new Worker(workerUrl,{type:'module',name:'matrix-geometry'});worker=current;
      get('geometry-compute').textContent='Cancel computation';get('geometry-play').disabled=true;get('geometry-progress').disabled=true;
      current.onmessage=({data:message})=>{
        if(worker!==current)return;
        if(message.type==='progress'){status(message.message);return;}
        if(message.type==='error'){status(`Geometry could not run: ${message.message}`,'failed');endWorker();return;}
        if(message.type!=='geometry')return;
        scene=message.result;progress=0;step=0;
        get('geometry-mode').value='morph';get('geometry-mode').options[1].disabled=!scene.trace_available;
        get('geometry-play').disabled=false;get('geometry-progress').disabled=false;get('geometry-reset').disabled=false;
        get('geometry-step').max=scene.steps.length;get('geometry-step-controls').hidden=true;get('geometry-progress').parentElement.hidden=false;
        get('geometry-badge').textContent='VERIFIED · WASM';
        status(`${scene.shared_checks.passed}/${scene.shared_checks.total} shared checks passed. Matrix and 125 field vectors independently verified.`,'pass');
        get('geometry-trace-note').textContent=scene.trace_available
          ? `${scene.steps.length} arithmetic updates captured inside the instrumented C kernel. Snapshots follow row, column, then inner-product index. Playback duration is illustrative.`
          : scene.operation==='multiply' ? (scene.trace_reason || bundle.trace_reason || 'This bundle has no step instrumentation. Build the optional trace WASM module to enable it.') : 'Calculation steps currently cover multiplication. This operation still uses verified WASM results.';
        get('geometry-profile-link').disabled=scene.operation!=='multiply'||!hasProfile();
        endWorker();render();
      };
      current.onerror=event=>{if(worker!==current)return;event.preventDefault();status('The browser could not initialize the geometry worker.','failed');endWorker();};
      timeout=setTimeout(()=>{if(worker===current){endWorker();status('Stopped after 30 seconds. Try computing again.','failed');}},30000);
      current.postMessage({type:'geometry',bundle,config});
    }catch(error){endWorker();status(`Geometry could not start: ${error.message}`,'failed');}
  }
  function hasProfile(){const data=JSON.parse(get('data').textContent);return data.profiles.some(p=>p.implementation==='wasm'&&p.status==='available'&&p.stacks.some(s=>s.some(f=>/^(?:wm_multiply|m_multiply)(?:[ (]|$)/.test(f))));}
  function operator(){if(!scene)return identity;if(isSteps())return step?scene.steps[step-1].matrix:Array(9).fill(0);if(scene.rotation_frames.length)return scene.rotation_frames[Math.round(progress*(scene.rotation_frames.length-1))];return scene.m.map((v,i)=>v*(1-progress)+scene.result[i]*progress);}
  function inspector(a){
    const frame=scene&&isSteps()&&step?scene.steps[step-1]:null;
    cells.forEach((cell,i)=>{cell.textContent=scene?fmt(a[i]):'—';cell.classList.toggle('active',!!frame&&i===frame.row*3+frame.col);});
    for(const [name,inputs]of Object.entries(fields))inputs.forEach((input,i)=>input.classList.toggle('active',!!frame&&i===(name==='m'?frame.row*3+frame.k:frame.k*3+frame.col)));
    get('geometry-equation').textContent=scene?scene.operation==='rotate'?'Mx → R(u, θ)Mx':scene.operation==='cross'?'Mx → (Mx) × u':`Mx → (${symbols[scene.operation]})x`:'Mx → (MN)x';
    get('geometry-description').textContent=scene?.operation==='rotate'?'Arrows show R(u, tθ)Mx at fixed grid points x. Angle playback preserves lengths and signed volume; all sampled operators are computed and verified in WebAssembly.':'Arrows show Ax at fixed grid points x. The wire cube and colored basis axes show the same operator acting on geometry.';
    get('geometry-matrix-title').textContent=isSteps()?(step===scene?.steps.length?'Computed result R':'Partial result R + active accumulator'):'Current operator A(t)';
    const volume=determinant(a);get('geometry-volume').textContent=scene?fmt(volume):'—';
    get('geometry-state').textContent=!scene?'Awaiting computation':isSteps()?(step===scene.steps.length?'Complete product':'Partial calculation'):Math.abs(volume)<1e-10?'Near-zero volume':volume<0?'Orientation reversed':'Orientation preserved';
    get('geometry-progress-label').textContent=isSteps()?`STEP ${step} / ${scene?.steps.length||0}`:`t = ${progress.toFixed(2)}`;
    get('geometry-progress').value=String(Math.round(progress*1000));get('geometry-step').value=String(step);
    get('geometry-step-label').textContent=`${step} / ${scene?.steps.length||0} multiply-adds`;
    get('geometry-prev').disabled=!scene||step===0;get('geometry-next').disabled=!scene||step===scene.steps.length;
    get('geometry-term').textContent=frame?`R[${frame.row+1},${frame.col+1}] += M[${frame.row+1},${frame.k+1}] × N[${frame.k+1},${frame.col+1}]\n${fmt(frame.left)} × ${fmt(frame.right)} = ${fmt(frame.product)}; accumulator = ${fmt(frame.sum)}\nC indices: row=${frame.row}, col=${frame.col}, k=${frame.k}`
      :scene&&isSteps()?'The result starts at zero. Step forward to accumulate its first dot product.':scene?`det(M) = ${fmt(scene.determinant_m)} · det(R) = ${fmt(scene.determinant_result)}\nArrows use one fixed display scale throughout the morph.`:'Compute a verified field to begin.';
    const eigen = scene?.eigen_m;
    const complexText=(re,im)=>im===0?fmt(re):`${fmt(re)} ${im<0?'−':'+'} ${fmt(Math.abs(im))}i`;
    get('geometry-eigen-values').textContent = eigen ? eigen.values.map((value,col)=>
      `λ${col+1} = ${complexText(value,eigen.imag_values[col])}; v${col+1} = (${[0,1,2].map(row=>complexText(eigen.vectors[row*3+col],eigen.imag_vectors[row*3+col])).join(', ')})`).join('\n')
      :'Compute an input matrix to inspect its real or complex eigenpairs.';
    const legend=[];
    if(eigen&&visibility.axes.checked&&eigen.imag_values.some(value=>value===0))legend.push('Purple arrows: real eigenvectors of M.');
    if(eigen&&visibility.ellipse.checked&&eigen.imag_values.some(value=>value>0))legend.push('Gold ellipse: complex invariant plane.');
    get('geometry-eigen-legend').textContent=legend.join(' ');
    get('geometry-eigen-legend').hidden=legend.length===0;
    if(get('geometry-frames').dataset.step!==String(frame?step:0)){
      get('geometry-frames').dataset.step=String(frame?step:0);get('geometry-frames').replaceChildren();
      for(const name of frame?.call_path||[]){const item=document.createElement('li');item.textContent=name;get('geometry-frames').append(item);}
      for(const [line,row]of sourceLines)row.classList.toggle('active',line===frame?.source_line);
      const active=frame&&sourceLines.get(frame.source_line);if(active){const pre=get('geometry-source');pre.scrollTop=Math.max(0,active.offsetTop-pre.offsetTop-pre.clientHeight/2);}
    }
  }
  function render(updateInspector=true){
    const showAxes=visibility.axes.checked,showBox=visibility.box.checked,showEllipse=visibility.ellipse.checked,showGrid=visibility.grid.checked;
    const bounds=canvas.getBoundingClientRect(),w=Math.max(1,bounds.width),h=Math.max(1,bounds.height),dpr=Math.min(devicePixelRatio||1,2);
    if(canvas.width!==Math.round(w*dpr)||canvas.height!==Math.round(h*dpr)){canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);}
    ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);
    const a=operator();if(updateInspector)inspector(a);
    const rgb=arrowColors[appearance.color.value]||[Math.round(100+145*progress),Math.round(207-28*progress),Math.round(189-71*progress)];
    const fieldColor=`rgb(${rgb.join(',')})`,fieldWidth=Number(appearance.width.value);
    get('geometry-field-swatch').style.background=fieldColor;
    get('geometry-arrow-note').textContent=reducedMotion.matches && appearance.motion.value!=='static'
      ? 'Your reduced-motion preference keeps arrow effects static. Play transformation and the progress sliders remain available.'
      : 'Effects decorate the fixed vectors, even while the matrix is paused. Faint teal arrows mark the input field; overlay colors stay fixed.';
    get('geometry-orbit-note').textContent=reducedMotion.matches && autoOrbit.checked
      ? 'Auto-orbit is paused by your reduced-motion preference. Drag or use arrow keys to orbit.'
      : '30 seconds per orbit. Drag to pause; release to resume.';
    const rowBounds=matrix=>[0,1,2].map(r=>Math.abs(matrix[r*3])+Math.abs(matrix[r*3+1])+Math.abs(matrix[r*3+2]));
    // Field arrows are normalized independently of the transformed overlays.
    // Fit visible overlays across the whole animation so playback keeps a steady scale.
    let extent=2;
    if(scene&&(showBox||showAxes)){
      extent=Math.max(extent,scene.rotation_extent||0);
      const matrices=[scene.m,scene.result,...(isSteps()?scene.steps.map(s=>s.matrix):[])];
      for(const matrix of matrices)extent=Math.max(extent,...rowBounds(matrix).map(x=>.55*x));
    }
    const fit=Math.min(w,h)*.34/extent*zoom;
    const project=v=>{const x=v[0]*Math.cos(yaw)+v[2]*Math.sin(yaw),z=-v[0]*Math.sin(yaw)+v[2]*Math.cos(yaw),y=v[1]*Math.cos(pitch)-z*Math.sin(pitch),depth=v[1]*Math.sin(pitch)+z*Math.cos(pitch),p=1/(1+depth/(extent*7));return [w/2+x*fit*p,h/2-y*fit*p,depth];};
    const lines=[];
    const line=(from,to,color,width=1,arrow=false,field=false)=>{const p=project(from),q=project(to);lines.push({p,q,color,width,arrow,field,depth:(p[2]+q[2])/2});};
    if(showGrid)for(let i=-2;i<=2;i+=.5){line([-2,-1.1,i],[2,-1.1,i],'rgba(121,166,179,.13)');line([i,-1.1,-2],[i,-1.1,2],'rgba(121,166,179,.13)');}
    const basis=[[1,0,0],[0,1,0],[0,0,1]],colors=['#fa8e88','#a0d790','#8ab8f8'];
    if(showAxes)basis.forEach(b=>line(scale(b,-1.55),scale(b,1.8),'rgba(179,204,210,.24)',1,true));
    if(scene){
      const maxNorm=Math.max(1,...scene.input_vectors.map(v=>Math.hypot(...v)),...scene.output_vectors.map(v=>Math.hypot(...v)),...(isSteps()?scene.steps.map(s=>Math.hypot(...rowBounds(s.matrix))):[]));
      const arrowScale=.72/maxNorm;
      scene.points.forEach((p,i)=>{
        const initial=scene.input_vectors[i],v=isSteps()||scene.rotation_frames.length?transform(a,p):initial.map((x,j)=>x*(1-progress)+scene.output_vectors[i][j]*progress);
        if(Math.hypot(...initial)>1e-10)line(p,add(p,scale(initial,arrowScale)),'rgba(100,207,189,.14)',.8,true);
        if(Math.hypot(...v)>1e-10)line(p,add(p,scale(v,arrowScale)),fieldColor,fieldWidth,true,true);
        else {const q=project(p);lines.push({p:q,q:[q[0]+1,q[1]+1,q[2]],color:fieldColor,width:fieldWidth+1,depth:q[2]});}
      });
      if(showBox){
        const cube=Array.from({length:8},(_,i)=>[(i&1)?.5:-.5,(i&2)?.5:-.5,(i&4)?.5:-.5]);
        cube.forEach((v,i)=>{for(const bit of [1,2,4])if(!(i&bit)){line(transform(scene.m,v),transform(scene.m,cube[i|bit]),'rgba(100,207,189,.23)');line(transform(a,v),transform(a,cube[i|bit]),'#f5b376',1.7);}});
      }
      if(showAxes)basis.forEach((b,i)=>line([0,0,0],transform(a,b),colors[i],2.7,true));
      if(scene.eigen_m) for(let col=0;col<3;col++) {
        const v=[0,1,2].map(row=>scene.eigen_m.vectors[row*3+col]);
        if(scene.eigen_m.imag_values[col]===0&&showAxes) line([0,0,0],v,'#d1a4ff',3,true);
        else if(scene.eigen_m.imag_values[col]>0&&showEllipse) {
          const imaginary=[0,1,2].map(row=>scene.eigen_m.imag_vectors[row*3+col]);
          const planePoint=angle=>v.map((x,i)=>1.6*(x*Math.cos(angle)+imaginary[i]*Math.sin(angle)));
          for(let segment=0;segment<64;segment++) line(planePoint(segment*Math.PI/32),planePoint((segment+1)*Math.PI/32),'#f1bd70',2);
          if(showAxes){
            line(scale(v,-1.4),scale(v,1.4),'#f1bd70',1.4);
            line(scale(imaginary,-1.4),scale(imaginary,1.4),'#f1bd70',1.4);
          }
        }
      }
    }
    lines.sort((a,b)=>b.depth-a.depth).forEach(({p,q,color,width,arrow,field})=>{
      const dx=q[0]-p[0],dy=q[1]-p[1],length=Math.hypot(dx,dy);
      const style=field?appearance.style.value:'open',motion=field&&arrowsMoving()?appearance.motion.value:'static';
      ctx.save();ctx.strokeStyle=color;ctx.fillStyle=color;ctx.lineWidth=width;ctx.lineCap='round';ctx.lineJoin='round';
      // Pulse changes opacity only. Every style keeps the same base and tip.
      if(field)ctx.globalAlpha=motion==='pulse'?.55+.4*(.5+.5*Math.sin(arrowPhase*Math.PI*2)):.85;
      const shaft=()=>{ctx.beginPath();ctx.moveTo(p[0],p[1]);ctx.lineTo(q[0],q[1]);ctx.stroke();};
      if(!arrow||length<=2)shaft();
      else {
        const n=Math.min(field?5+width:5,length*.35),ux=dx/length,uy=dy/length;
        const at=(back,side)=>[q[0]-back*ux-side*uy,q[1]-back*uy+side*ux];
        if(style==='tapered'){
          ctx.beginPath();ctx.moveTo(p[0],p[1]);
          ctx.lineTo(...at(n,n*.22));ctx.lineTo(...at(n,n*.6));ctx.lineTo(q[0],q[1]);
          ctx.lineTo(...at(n,-n*.6));ctx.lineTo(...at(n,-n*.22));ctx.closePath();ctx.fill();
        }else{
          shaft();ctx.beginPath();ctx.moveTo(...at(n,n*.5));ctx.lineTo(q[0],q[1]);ctx.lineTo(...at(n,-n*.5));
          if(style==='solid'){ctx.closePath();ctx.fill();}else ctx.stroke();
          if(style==='chevron'){
            ctx.beginPath();ctx.moveTo(...at(n*1.8,n*.5));ctx.lineTo(...at(n*.8,0));ctx.lineTo(...at(n*1.8,-n*.5));ctx.stroke();
          }
        }
        if(motion==='flow'){
          // Travel along the shaft without moving the vector endpoint or length.
          const radius=Math.min(width+.65,length*.08),distance=radius+(length-2*radius)*arrowPhase;
          ctx.globalAlpha=.95;ctx.fillStyle='#ffffff';ctx.beginPath();
          ctx.arc(p[0]+ux*distance,p[1]+uy*distance,radius,0,Math.PI*2);ctx.fill();
        }
      }
      ctx.restore();
    });
    if(showAxes){ctx.font='11px ui-monospace,monospace';basis.forEach((b,i)=>{const p=project(scale(b,1.94));ctx.fillStyle=colors[i];ctx.fillText(['x','y','z'][i],p[0],p[1]);});}
    if(!scene){ctx.fillStyle='#a8c1c9';ctx.font='13px system-ui';ctx.textAlign='center';ctx.fillText('Compute a matrix field to explore it here',w/2,h*.82);ctx.textAlign='left';}
    syncAnimation();
  }
  function tick(time){
    animation=null;
    if(!scene||document.hidden){previousTime=null;arrowPreviousTime=null;orbitPreviousTime=null;return;}
    const transforming=playing;
    if(playing){
      if(previousTime!==null){const delta=Math.min(100,time-previousTime)/(Number(get('geometry-duration').value)*1000);progress=Math.min(1,progress+delta);if(isSteps())step=Math.min(scene.steps.length,Math.floor(progress*scene.steps.length));}
      previousTime=time;if(progress>=1)stopAnimation();
    }
    if(arrowsMoving()){
      if(arrowPreviousTime!==null)arrowPhase=(arrowPhase+Math.min(100,time-arrowPreviousTime)/1800)%1;
      arrowPreviousTime=time;
    }else arrowPreviousTime=null;
    if(orbitMoving()){
      if(orbitPreviousTime!==null)yaw=(yaw+Math.min(100,time-orbitPreviousTime)*Math.PI*2/30000)%(Math.PI*2);
      orbitPreviousTime=time;
    }else orbitPreviousTime=null;
    render(transforming);
  }
  get('geometry-play').onclick=()=>{if(!scene)return;if(playing){stopAnimation();return;}if(progress>=1){progress=0;step=0;}playing=true;previousTime=null;get('geometry-play').textContent='Pause animation';syncAnimation();};
  get('geometry-progress').oninput=()=>{stopAnimation();progress=Number(get('geometry-progress').value)/1000;render();};
  function selectStep(value){stopAnimation();step=Math.max(0,Math.min(scene?.steps.length||0,value));progress=scene?.steps.length?step/scene.steps.length:0;render();}
  get('geometry-step').oninput=()=>selectStep(Number(get('geometry-step').value));
  get('geometry-prev').onclick=()=>selectStep(step-1);get('geometry-next').onclick=()=>selectStep(step+1);
  get('geometry-mode').onchange=()=>{stopAnimation();progress=0;step=0;get('geometry-step-controls').hidden=!isSteps();get('geometry-progress').parentElement.hidden=isSteps();if(isSteps())get('geometry-trace-details').open=true;get('geometry-play').textContent=isSteps()?'Play calculation':'Play transformation';render();};
  get('geometry-reset').onclick=()=>{stopAnimation();progress=0;step=0;yaw=-.68;pitch=.43;zoom=1;orbitPreviousTime=null;render();};
  get('geometry-preset').onchange=()=>{applyPreset();if(get('geometry-preset').value!=='custom')compute();};
  get('geometry-randomize').onclick=()=>{
    for(const inputs of Object.values(fields))for(const input of inputs)input.value=(Math.floor(Math.random()*17)-8)/4;
    get('geometry-preset').value='custom';
    invalidate();
    if(!get('geometry-compute').disabled)compute();
    else status('Matrices randomized. WebAssembly is unavailable in this report.','unavailable');
  };
  get('geometry-operation').onchange=operationChanged;get('geometry-scalar').oninput=invalidate;
  for(const input of [...axisInputs,get('geometry-angle')])input.oninput=()=>{get('geometry-preset').value='custom';invalidate();};
  get('geometry-compute').onclick=compute;
  for(const input of Object.values(visibility))input.onchange=()=>render();
  for(const input of Object.values(appearance))input.onchange=()=>{arrowPreviousTime=null;render();};
  autoOrbit.onchange=()=>{orbitPreviousTime=null;render();};
  get('geometry-profile-link').onclick=()=>{stopAnimation();window.dispatchEvent(new Event('matrix-profile-focus'));};
  canvas.onpointerdown=e=>{drag=[e.clientX,e.clientY];orbitPreviousTime=null;canvas.setPointerCapture(e.pointerId);syncAnimation();};
  canvas.onpointermove=e=>{if(!drag)return;yaw+=(e.clientX-drag[0])*.009;pitch=Math.max(-1.3,Math.min(1.3,pitch+(e.clientY-drag[1])*.009));drag=[e.clientX,e.clientY];render();};
  canvas.onpointerup=canvas.onpointercancel=canvas.onlostpointercapture=()=>{drag=null;orbitPreviousTime=null;syncAnimation();};
  canvas.addEventListener('wheel',e=>{e.preventDefault();changeZoom(Math.exp(-e.deltaY*.001));},{passive:false});
  canvas.onkeydown=e=>{if(!['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','+','-'].includes(e.key))return;e.preventDefault();if(e.key==='+'){changeZoom(1.1);return;}if(e.key==='-'){changeZoom(1/1.1);return;}if(e.key==='ArrowLeft')yaw-=.12;if(e.key==='ArrowRight')yaw+=.12;if(e.key==='ArrowUp')pitch=Math.min(1.3,pitch+.12);if(e.key==='ArrowDown')pitch=Math.max(-1.3,pitch-.12);render();};
  document.addEventListener('visibilitychange',()=>{arrowPreviousTime=null;orbitPreviousTime=null;if(document.hidden)stopAnimation();else render();});
  reducedMotion.addEventListener('change',()=>{arrowPreviousTime=null;orbitPreviousTime=null;render();});
  new ResizeObserver(()=>render()).observe(canvas);
  applyPreset();
  if(!bundle.available||!bundle.geometry_worker_source||typeof Worker==='undefined'||typeof WebAssembly==='undefined')status(bundle.reason||'This report has no geometry worker. Regenerate it with a current WASM build.','unavailable');
  else {get('geometry-compute').disabled=false;compute();}
})();
