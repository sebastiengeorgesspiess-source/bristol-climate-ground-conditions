(function(root){
'use strict';
const mean=a=>a.reduce((s,x)=>s+x,0)/a.length;
const sd=a=>Math.sqrt(a.reduce((s,x)=>s+(x-mean(a))**2,0)/(a.length-1));
const month=d=>+d.slice(0,4)*12 + +d.slice(5,7)-1;
function corr(a,b){if(a.length<3||a.length!==b.length)return null;const ma=mean(a),mb=mean(b);let n=0,da=0,db=0;for(let i=0;i<a.length;i++){n+=(a[i]-ma)*(b[i]-mb);da+=(a[i]-ma)**2;db+=(b[i]-mb)**2;}return da>0&&db>0?n/Math.sqrt(da*db):null;}
function rank(a){const order=a.map((v,i)=>({v,i})).sort((a,b)=>a.v-b.v),out=[];for(let i=0;i<order.length;){let j=i+1;while(j<order.length&&order[j].v===order[i].v)j++;for(let k=i;k<j;k++)out[order[k].i]=(i+1+j)/2;i=j;}return out;}
function transform(rows,mode){if(mode==='raw')return rows;const base={};rows.filter(r=>r.date>='1991-01'&&r.date<='2020-12').forEach(r=>(base[r.date.slice(5)]??=[]).push(r.value));return rows.flatMap(r=>{const b=base[r.date.slice(5)]||[],s=sd(b);if(b.length<24||(mode==='z'&&!(s>0)))return [];return [{date:r.date,value:(r.value-mean(b))/(mode==='z'?s:1)}];});}
function smooth(rows,n){return rows.flatMap((r,i)=>{const w=rows.slice(i-n+1,i+1);return i>=n-1&&w.every((x,j)=>month(x.date)===month(r.date)-n+1+j)?[{date:r.date,value:mean(w.map(x=>x.value))}]:[];});}
function pairs(a,b,mode,n,lag,from,to){const bm=new Map(smooth(transform(b,mode),n).map(x=>[month(x.date),x]));return smooth(transform(a,mode),n).flatMap(x=>{const y=bm.get(month(x.date)+lag);return y&&x.date>=from+'-01'&&x.date<=to+'-12'&&y.date>=from+'-01'&&y.date<=to+'-12'?[{date:x.date,dateB:y.date,a:x.value,b:y.value}]:[];});}
function cdf(x){const t=1/(1+.2316419*Math.abs(x)),d=.3989423*Math.exp(-x*x/2),p=1-d*t*(.31938153+t*(-.356563782+t*(1.781477937+t*(-1.821255978+t*1.330274429))));return x>=0?p:1-p;}
function infer(rows){const a=rows.map(x=>x.a),b=rows.map(x=>x.b),r=corr(a,b);if(r===null||rows.length<24)return null;const adj=rows.slice(1).map((x,i)=>({x,prev:rows[i]})).filter(z=>month(z.x.date)-month(z.prev.date)===1);if(adj.length<12)return null;const ra=corr(adj.map(z=>z.prev.a),adj.map(z=>z.x.a)),rb=corr(adj.map(z=>z.prev.b),adj.map(z=>z.x.b));if(ra===null||rb===null)return null;const ne=Math.min(rows.length,rows.length*(1-ra*rb)/(1+ra*rb));if(ne<=3)return null;const z=Math.atanh(Math.max(-.999999,Math.min(.999999,r))),se=1/Math.sqrt(ne-3);return {r,n:rows.length,ne,lo:Math.tanh(z-1.96*se),hi:Math.tanh(z+1.96*se),p:Math.max(0,Math.min(1,2*(1-cdf(Math.abs(z)/se))))};}
function bh(tests){const sorted=tests.map((t,i)=>({...t,i})).sort((a,b)=>a.p-b.p);let q=1;const out=[];for(let i=sorted.length-1;i>=0;i--){q=Math.min(q,sorted[i].p*sorted.length/(i+1));out[sorted[i].i]=q;}return out;}
const api={mean,sd,month,corr,rank,transform,smooth,pairs,infer,bh};if(typeof module!=='undefined')module.exports=api;else root.BristolMath=api;
})(typeof window!=='undefined'?window:this);
