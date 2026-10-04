const $=id=>document.getElementById(id),bundle=JSON.parse($('bundle').textContent),states={},colors=()=>{const s=getComputedStyle(document.documentElement);return ['--ink','--muted','--green','--orange','--plot'].map(v=>s.getPropertyValue(v).trim());};
const number=v=>v===null?'undefined':Number.isFinite(v)?v.toPrecision(5):'outside range';
$('provenance').textContent=`Embedded bundle SHA-256: ${bundle.sha256}`;
addEventListener('matrix-theme-change',()=>{for(const kind of ['pca','fit'])if(states[kind])draw(kind);});
function stop(s){s.worker?.terminate();s.worker=null;if(s.url)URL.revokeObjectURL(s.url);s.url=null;clearTimeout(s.watchdog);}
function imageOutputState(state){const s=states.image,panel=$('image-output-panel');panel.dataset.state=state;panel.setAttribute('aria-busy',String(state==='updating'));$('image-output-state').textContent=state==='ready'?`Showing rank ${s.displayedRank}.`:state==='updating'?(s.displayedRank===null?'Updating reconstruction…':`Updating to rank ${$('image-rank').value}; showing previous rank ${s.displayedRank}.`):(s.displayedRank===null?'No reconstruction available.':`Update failed; showing previous rank ${s.displayedRank}.`);}
function failure(kind,message){const s=states[kind];++s.id;clearTimeout(s.watchdog);s.result=null;$(`${kind}-status`).textContent=message;if(kind==='image')imageOutputState('error');}
function schedule(kind,{rankOnly=false}={}){const s=states[kind];++s.id;clearTimeout(s.timer);clearTimeout(s.watchdog);s.result=null;$(`${kind}-metrics`).textContent='';$(`${kind}-status`).textContent='Waiting for editing to pause…';if(kind!=='image')draw(kind);else{
 if(!rankOnly){s.displayedRank=null;s.spectrumValues=null;$('image-output').getContext('2d').clearRect(0,0,24,24);$('image-spectrum').replaceChildren();$('image-spectrum-values').replaceChildren();$('image-spectrum-cutoff').hidden=true;}
 if(s.spectrumValues)spectrum(s.spectrumValues);else $('image-spectrum-caption').textContent=`Rank ${$('image-rank').value} selected. Calculating singular values…`;
 imageOutputState('updating');
 }
 s.timer=setTimeout(()=>compute(kind),300);
}
function compute(kind){const s=states[kind],id=s.id;
 try{if(kind!=='image'&&['x','y'].some(axis=>{const raw=$(`${kind}-${axis}`).value;return raw.trim()===''||!Number.isFinite(Number(raw))||Math.abs(Number(raw))>5;}))throw new RangeError('Enter coordinates between −5 and 5.');if(!s.worker){s.url=URL.createObjectURL(new Blob([bundle.worker_source],{type:'text/javascript'}));s.worker=new Worker(s.url,{type:'module'});s.worker.onmessage=({data})=>{if(data.id!==s.id)return;clearTimeout(s.watchdog);if(data.error){failure(kind,data.error);return;}s.result=data.result;show(kind);};s.worker.onerror=()=>{stop(s);failure(kind,'Calculation could not run. Change an input to retry.');};}
 $(`${kind}-status`).textContent='Calculating in WebAssembly…';s.watchdog=setTimeout(()=>{stop(s);failure(kind,'Calculation timed out. Change an input to retry.');},20000);
 s.worker.postMessage({id,kind,bundle,config:kind==='image'?{pixels:s.pixels,size:24,rank:Number($('image-rank').value)}:{points:s.points,degree:Number($('fit-degree').value)}});
 }catch(error){stop(s);failure(kind,error.message);}
}
function show(kind){const s=states[kind],r=s.result;
 if(kind==='pca'){$('pca-metrics').textContent=`Mean: [${r.mean.map(number).join(', ')}]\nVariances: ${r.values.map(number).join(', ')}\nFirst component: ${r.explained===null?'undefined':(100*r.explained).toFixed(2)+'%'} of variance`;$('pca-status').textContent=r.ambiguous?'No unique principal direction: variance is zero or tied.':`${s.points.length} points · sample covariance → symmetric eigenvectors`;draw(kind);}
 if(kind==='fit'){$('fit-metrics').textContent=`Coefficients for powers of t=x/5:\n${r.coefficients.map(number).join(', ')}\nRMSE: ${number(r.rmse)}\nR²: ${number(r.r2)}`;$('fit-status').textContent=`${s.points.length} observations · column-pivoted QR`;draw(kind);}
 if(kind==='image'){paint($('image-output'),r.reconstructed);s.displayedRank=Number($('image-rank').value);s.spectrumValues=r.singularValues;imageOutputState('ready');$('image-metrics').textContent=`Relative error: ${(100*r.relativeError).toFixed(2)}%\nRMSE: ${number(r.rmse)}\nRetained energy: ${r.retained===null?'undefined':(100*r.retained).toFixed(2)+'%'}\nScalars: ${r.factorScalars} / ${r.originalScalars}\nFactor storage: ${(100*r.factorScalars/r.originalScalars).toFixed(1)}% of original`;$('image-status').textContent='24 × 24 grayscale · WASM one-sided Jacobi SVD';spectrum(r.singularValues);}
}
function spectrum(values){
 const rank=Number($('image-rank').value),max=Math.max(...values,1e-30),total=values.length;
 $('image-spectrum').replaceChildren(...values.map((value,i)=>{const bar=document.createElement('span');bar.style.height=`${100*value/max}%`;bar.className=i<rank?'retained':'omitted';bar.title=`σ${i+1} = ${number(value)} · ${i<rank?'retained':'omitted'}`;return bar;}));
 const cutoff=$('image-spectrum-cutoff');cutoff.hidden=!total;cutoff.className='spectrum-cutoff '+(rank===0?'at-start':rank===total?'at-end':'');cutoff.style.left=rank===0?'0%':rank===total?'100%':`calc(${100*rank/total}% + ${3*(rank/total-.5)}px)`;
 $('image-spectrum-cutoff-label').textContent=`k = ${rank}`;
 $('image-spectrum-caption').textContent=rank===0?`Selected rank 0: all ${total} components omitted.`:rank===total?`Selected rank ${rank}: all ${total} components retained.`:`Selected rank ${rank}: retain σ1–σ${rank}; omit σ${rank+1}–σ${total}.`;
 $('image-spectrum-values').replaceChildren(...values.map((value,i)=>{const row=document.createElement('tr'),name=document.createElement('th'),cell=document.createElement('td'),retained=document.createElement('td'),included=i<rank;name.scope='row';name.textContent=`σ${i+1}`;cell.textContent=number(value);retained.textContent=included?'Retained':'Omitted';row.className=included?'retained':'';row.append(name,cell,retained);return row;}));
}
function draw(kind){
 const s=states[kind],canvas=$(`${kind}-canvas`),ctx=canvas.getContext('2d'),w=canvas.clientWidth,h=canvas.clientHeight,dpr=Math.min(2,devicePixelRatio||1),[ink,muted,green,orange,plot]=colors();
 canvas.width=w*dpr;canvas.height=h*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);
 const unit=Math.min(w,h)/10,xy=([x,y])=>[w/2+x*unit,h/2-y*unit],labels=[];
 ctx.strokeStyle=muted;ctx.lineWidth=.5;ctx.globalAlpha=.3;
 for(let i=-4;i<=4;i++){if(i===0)continue;ctx.beginPath();ctx.moveTo(...xy([i,-5]));ctx.lineTo(...xy([i,5]));ctx.moveTo(...xy([-5,i]));ctx.lineTo(...xy([5,i]));ctx.stroke();}
 ctx.globalAlpha=1;ctx.strokeStyle=muted;ctx.lineWidth=1.5;ctx.beginPath();ctx.moveTo(...xy([-5,0]));ctx.lineTo(...xy([5,0]));ctx.moveTo(...xy([0,-5]));ctx.lineTo(...xy([0,5]));ctx.stroke();
 ctx.fillStyle=muted;ctx.font='13px system-ui';
 for(const value of [-4,-2,2,4]){const [x,y]=xy([value,value]);ctx.beginPath();ctx.moveTo(x,h/2-3);ctx.lineTo(x,h/2+3);ctx.moveTo(w/2-3,y);ctx.lineTo(w/2+3,y);ctx.stroke();ctx.textAlign='center';ctx.fillText(String(value),x,h/2+18);ctx.textAlign='right';ctx.fillText(String(value),w/2-9,y+4);}
 ctx.textAlign='left';ctx.fillText('0',w/2+7,h/2+18);ctx.fillStyle=ink;ctx.font='600 14px system-ui';ctx.fillText('x',Math.min(w-16,w/2+5*unit-12),h/2-10);ctx.fillText('y',w/2+10,Math.max(18,h/2-5*unit+16));
 const line=(points,color,dash=[])=>{ctx.strokeStyle=color;ctx.lineWidth=2;ctx.setLineDash(dash);ctx.beginPath();points.forEach((p,i)=>i?ctx.lineTo(...xy(p)):ctx.moveTo(...xy(p)));ctx.stroke();ctx.setLineDash([]);};
 if(s.result&&kind==='pca'){
  const r=s.result;r.projected.forEach((p,i)=>line([s.points[i],p],muted,[2,4]));
  if(!r.ambiguous)r.axes.forEach((axis,j)=>{const length=2*Math.sqrt(r.values[j]),endpoints=[-1,1].map(sign=>r.mean.map((v,i)=>v+sign*axis[i]*length));line(endpoints,j?orange:green);const end=axis[0]<0||(axis[0]===0&&axis[1]<0)?endpoints[0]:endpoints[1];labels.push({text:`PC${j+1}`,point:xy(end),origin:xy(r.mean),color:j?orange:green});});
 }
 if(s.result&&kind==='fit'){const r=s.result;r.predicted.forEach((y,i)=>line([s.points[i],[s.points[i][0],y]],orange,[3,3]));line(r.curve,green);}
 s.points.forEach((p,i)=>{ctx.beginPath();ctx.arc(...xy(p),i===s.selected?6:4,0,Math.PI*2);ctx.fillStyle=i===s.selected?orange:ink;ctx.fill();});
 // Anchor tags to the visible end of each direction, including clipped axes.
 const placed=[],clamp=(v,lo,hi)=>Math.max(lo,Math.min(hi,v));ctx.font='600 13px system-ui';ctx.textAlign='left';
 for(const {text,point,origin,color} of labels){
  const delta=point.map((v,i)=>v-origin[i]);let t=1;
  for(let i=0;i<2;i++){const extent=i?h:w;if(point[i]<8&&delta[i]<0)t=Math.min(t,(8-origin[i])/delta[i]);if(point[i]>extent-8&&delta[i]>0)t=Math.min(t,(extent-8-origin[i])/delta[i]);}
  const anchor=origin.map((v,i)=>clamp(v+delta[i]*clamp(t,0,1),8,(i?h:w)-8)),width=ctx.measureText(text).width+12,height=24;
  let x=clamp(anchor[0]+8,6,w-width-6),y=clamp(anchor[1]-height-7,6,h-height-6);
  for(const rect of placed)if(x<rect.x+rect.width+4&&x+width+4>rect.x&&y<rect.y+height+4&&y+height+4>rect.y)y=clamp(rect.y+height+5<=h-height-6?rect.y+height+5:rect.y-height-5,6,h-height-6);
  ctx.strokeStyle=color;ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(...anchor);ctx.lineTo(clamp(anchor[0],x,x+width),clamp(anchor[1],y,y+height));ctx.stroke();ctx.fillStyle=plot;ctx.fillRect(x,y,width,height);ctx.fillStyle=color;ctx.fillText(text,x+6,y+17);placed.push({x,y,width});
 }
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
states.image={id:0,pixels:[],displayedRank:null,spectrumValues:null};
function imagePreset(){states.image.pixels=Array.from({length:576},(_,i)=>{const x=i%24,y=Math.floor(i/24),p=$('image-preset').value;return p==='blank'?0:p==='waves'?.5+.25*Math.sin(x/3)+.25*Math.cos(y/4):p==='checker'?((Math.floor(x/3)+Math.floor(y/3))%2):Math.hypot(x-8,y-9)<6?.9:x>12&&y>12?.65:.08;});paint($('image-original'),states.image.pixels);schedule('image');}
$('image-preset').onchange=imagePreset;$('image-rank').oninput=()=>{$('image-rank-label').textContent=$('image-rank').value;schedule('image',{rankOnly:true});};
let painting=false;const original=$('image-original');function brush(e){const r=original.getBoundingClientRect(),x=Math.max(0,Math.min(23,Math.floor((e.clientX-r.left)/r.width*24))),y=Math.max(0,Math.min(23,Math.floor((e.clientY-r.top)/r.height*24)));states.image.pixels[y*24+x]=Number($('image-ink').value);paint(original,states.image.pixels);schedule('image');}
original.onpointerdown=e=>{painting=true;original.setPointerCapture(e.pointerId);brush(e);};original.onpointermove=e=>{if(painting)brush(e);};original.onpointerup=original.onpointercancel=()=>painting=false;
imagePreset();addEventListener('resize',()=>{draw('pca');draw('fit');});addEventListener('pagehide',()=>{for(const s of Object.values(states)){clearTimeout(s.timer);stop(s);}});
