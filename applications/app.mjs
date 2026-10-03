const $=id=>document.getElementById(id),bundle=JSON.parse($('bundle').textContent),states={},colors=()=>{const s=getComputedStyle(document.documentElement);return ['--ink','--muted','--green','--orange','--plot'].map(v=>s.getPropertyValue(v).trim());};
const number=v=>v===null?'undefined':Number.isFinite(v)?v.toPrecision(5):'outside range';
$('provenance').textContent=`Embedded bundle SHA-256: ${bundle.sha256}`;
const media=matchMedia('(prefers-color-scheme: dark)');
try{$('theme').value=localStorage.getItem('linear-a-report-theme')||'system';}catch{}
function theme(){document.documentElement.dataset.theme=$('theme').value==='system'?(media.matches?'dark':'light'):$('theme').value;try{localStorage.setItem('linear-a-report-theme',$('theme').value);}catch{};for(const kind of ['pca','fit'])if(states[kind])draw(kind);}
$('theme').onchange=theme;media.addEventListener('change',theme);theme();
function stop(s){s.worker?.terminate();s.worker=null;if(s.url)URL.revokeObjectURL(s.url);s.url=null;clearTimeout(s.watchdog);}
function schedule(kind){const s=states[kind];++s.id;clearTimeout(s.timer);clearTimeout(s.watchdog);s.result=null;$(`${kind}-metrics`).textContent='';$(`${kind}-status`).textContent='Waiting for editing to pause…';if(kind!=='image')draw(kind);else{$('image-output').getContext('2d').clearRect(0,0,24,24);$('image-spectrum').replaceChildren();}
 s.timer=setTimeout(()=>compute(kind),300);
}
function compute(kind){const s=states[kind],id=s.id;
 try{if(kind!=='image'&&['x','y'].some(axis=>{const raw=$(`${kind}-${axis}`).value;return raw.trim()===''||!Number.isFinite(Number(raw))||Math.abs(Number(raw))>5;}))throw new RangeError('Enter coordinates between −5 and 5.');if(!s.worker){s.url=URL.createObjectURL(new Blob([bundle.worker_source],{type:'text/javascript'}));s.worker=new Worker(s.url,{type:'module'});s.worker.onmessage=({data})=>{if(data.id!==s.id)return;clearTimeout(s.watchdog);if(data.error){s.result=null;$(`${kind}-status`).textContent=data.error;return;}s.result=data.result;show(kind);};s.worker.onerror=()=>{stop(s);$(`${kind}-status`).textContent='Calculation could not run. Change an input to retry.';};}
 $(`${kind}-status`).textContent='Calculating in WebAssembly…';s.watchdog=setTimeout(()=>{++s.id;stop(s);$(`${kind}-status`).textContent='Calculation timed out. Change an input to retry.';},20000);
 s.worker.postMessage({id,kind,bundle,config:kind==='image'?{pixels:s.pixels,size:24,rank:Number($('image-rank').value)}:{points:s.points,degree:Number($('fit-degree').value)}});
 }catch(error){stop(s);$(`${kind}-status`).textContent=error.message;}
}
function show(kind){const s=states[kind],r=s.result;
 if(kind==='pca'){$('pca-metrics').textContent=`Mean: [${r.mean.map(number).join(', ')}]\nVariances: ${r.values.map(number).join(', ')}\nFirst component: ${r.explained===null?'undefined':(100*r.explained).toFixed(2)+'%'} of variance`;$('pca-status').textContent=r.ambiguous?'No unique principal direction: variance is zero or tied.':`${s.points.length} points · sample covariance → symmetric eigenvectors`;draw(kind);}
 if(kind==='fit'){$('fit-metrics').textContent=`Coefficients for powers of t=x/5:\n${r.coefficients.map(number).join(', ')}\nRMSE: ${number(r.rmse)}\nR²: ${number(r.r2)}`;$('fit-status').textContent=`${s.points.length} observations · column-pivoted QR`;draw(kind);}
 if(kind==='image'){paint($('image-output'),r.reconstructed);$('image-metrics').textContent=`Relative error: ${(100*r.relativeError).toFixed(2)}%\nRMSE: ${number(r.rmse)}\nRetained energy: ${r.retained===null?'undefined':(100*r.retained).toFixed(2)+'%'}\nScalars: ${r.factorScalars} / ${r.originalScalars}\nFactor storage: ${(100*r.factorScalars/r.originalScalars).toFixed(1)}% of original`;$('image-status').textContent='24 × 24 grayscale · WASM one-sided Jacobi SVD';const max=Math.max(...r.singularValues,1e-30);$('image-spectrum').replaceChildren(...r.singularValues.map((v,i)=>{const bar=document.createElement('span');bar.style.height=`${100*v/max}%`;bar.title=`σ${i+1} = ${number(v)}`;return bar;}));}
}
function draw(kind){const s=states[kind],canvas=$(`${kind}-canvas`),ctx=canvas.getContext('2d'),w=canvas.clientWidth,h=canvas.clientHeight,dpr=Math.min(2,devicePixelRatio||1),[ink,muted,green,orange]=colors();canvas.width=w*dpr;canvas.height=h*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);const unit=Math.min(w,h)/10,xy=([x,y])=>[w/2+x*unit,h/2-y*unit];
 ctx.strokeStyle=muted;ctx.lineWidth=.5;for(let i=-4;i<=4;i++){ctx.beginPath();ctx.moveTo(...xy([i,-5]));ctx.lineTo(...xy([i,5]));ctx.moveTo(...xy([-5,i]));ctx.lineTo(...xy([5,i]));ctx.stroke();}ctx.fillStyle=muted;ctx.font='11px system-ui';ctx.fillText('x / y · −5 … 5',10,16);
 const line=(points,color,dash=[])=>{ctx.strokeStyle=color;ctx.lineWidth=2;ctx.setLineDash(dash);ctx.beginPath();points.forEach((p,i)=>i?ctx.lineTo(...xy(p)):ctx.moveTo(...xy(p)));ctx.stroke();ctx.setLineDash([]);};
 if(s.result&&kind==='pca'){const r=s.result;r.projected.forEach((p,i)=>line([s.points[i],p],muted,[2,4]));if(!r.ambiguous)r.axes.forEach((a,j)=>{const length=2*Math.sqrt(r.values[j]);line([-1,1].map(sign=>r.mean.map((v,i)=>v+sign*a[i]*length)),j?orange:green);});}
 if(s.result&&kind==='fit'){const r=s.result;r.predicted.forEach((y,i)=>line([s.points[i],[s.points[i][0],y]],orange,[3,3]));line(r.curve,green);}
 s.points.forEach((p,i)=>{ctx.beginPath();ctx.arc(...xy(p),i===s.selected?6:4,0,Math.PI*2);ctx.fillStyle=i===s.selected?orange:ink;ctx.fill();});
}
function editor(kind){const s=states[kind],select=$(`${kind}-selected`);select.replaceChildren(...s.points.map((_,i)=>{const o=document.createElement('option');o.value=i;o.textContent=`Point ${i+1}`;return o;}));select.value=s.selected;$(`${kind}-x`).value=s.points[s.selected][0];$(`${kind}-y`).value=s.points[s.selected][1];}
for(const kind of ['pca','fit']){
 const preset=()=>Array.from({length:12},(_,i)=>{const x=-4+i*8/11;return [x,kind==='pca'?.6*x+Math.sin(i*2)*.75:.18*x*x-1+Math.sin(i*3)*.5];});
 const s=states[kind]={points:preset(),selected:0,id:0};editor(kind);
 $(`${kind}-selected`).onchange=()=>{s.selected=Number($(`${kind}-selected`).value);editor(kind);draw(kind);};
 for(const [axis,j] of [['x',0],['y',1]])$(`${kind}-${axis}`).oninput=()=>{const input=$(`${kind}-${axis}`),v=input.value.trim()===''?NaN:Number(input.value);if(!Number.isFinite(v)||Math.abs(v)>5){++s.id;clearTimeout(s.timer);clearTimeout(s.watchdog);s.result=null;$(`${kind}-metrics`).textContent='';$(`${kind}-status`).textContent='Enter coordinates between −5 and 5.';draw(kind);return;}s.points[s.selected][j]=v;schedule(kind);};
 $(`${kind}-reset`).onclick=()=>{s.points=preset();s.selected=0;editor(kind);schedule(kind);};
 $(`${kind}-remove`).onclick=()=>{if(s.points.length<=2)return;s.points.splice(s.selected,1);s.selected=Math.min(s.selected,s.points.length-1);editor(kind);schedule(kind);};
 $(`${kind}-add`).onclick=()=>{if(s.points.length>=80)return;s.points.push([0,0]);s.selected=s.points.length-1;editor(kind);schedule(kind);};
 const canvas=$(`${kind}-canvas`);let dragging=false;
 const locate=e=>{const r=canvas.getBoundingClientRect();const unit=Math.min(r.width,r.height)/10;return [Math.max(-5,Math.min(5,(e.clientX-r.left-r.width/2)/unit)),Math.max(-5,Math.min(5,(r.height/2-e.clientY+r.top)/unit))];};
 canvas.onpointerdown=e=>{const p=locate(e),distance=s.points.map(v=>Math.hypot(v[0]-p[0],v[1]-p[1])),nearest=distance.indexOf(Math.min(...distance));if(distance[nearest]<.35)s.selected=nearest;else if(s.points.length<80){s.points.push(p);s.selected=s.points.length-1;}else return;dragging=true;canvas.setPointerCapture(e.pointerId);editor(kind);schedule(kind);};
 canvas.onpointermove=e=>{if(!dragging)return;s.points[s.selected]=locate(e);editor(kind);schedule(kind);};canvas.onpointerup=canvas.onpointercancel=()=>dragging=false;schedule(kind);
}
$('fit-degree').onchange=()=>schedule('fit');
function paint(canvas,pixels){const ctx=canvas.getContext('2d'),image=ctx.createImageData(24,24);pixels.forEach((v,i)=>{const gray=Math.round(255*Math.max(0,Math.min(1,v)));image.data.set([gray,gray,gray,255],4*i);});ctx.putImageData(image,0,0);}
states.image={id:0,pixels:[]};
function imagePreset(){states.image.pixels=Array.from({length:576},(_,i)=>{const x=i%24,y=Math.floor(i/24),p=$('image-preset').value;return p==='blank'?0:p==='waves'?.5+.25*Math.sin(x/3)+.25*Math.cos(y/4):p==='checker'?((Math.floor(x/3)+Math.floor(y/3))%2):Math.hypot(x-8,y-9)<6?.9:x>12&&y>12?.65:.08;});paint($('image-original'),states.image.pixels);schedule('image');}
$('image-preset').onchange=imagePreset;$('image-rank').oninput=()=>{$('image-rank-label').textContent=$('image-rank').value;schedule('image');};
let painting=false;const original=$('image-original');function brush(e){const r=original.getBoundingClientRect(),x=Math.max(0,Math.min(23,Math.floor((e.clientX-r.left)/r.width*24))),y=Math.max(0,Math.min(23,Math.floor((e.clientY-r.top)/r.height*24)));states.image.pixels[y*24+x]=Number($('image-ink').value);paint(original,states.image.pixels);schedule('image');}
original.onpointerdown=e=>{painting=true;original.setPointerCapture(e.pointerId);brush(e);};original.onpointermove=e=>{if(painting)brush(e);};original.onpointerup=original.onpointercancel=()=>painting=false;
imagePreset();addEventListener('resize',()=>{draw('pca');draw('fit');});addEventListener('pagehide',()=>{for(const s of Object.values(states)){clearTimeout(s.timer);stop(s);}});
