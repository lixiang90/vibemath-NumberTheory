// Node.js standard library. Exact-arithmetic certificate for the fixed-core reductions.
function stage(run,input){
 let result;
 run(input,{log(message){
   if(result!==undefined)throw Error("duplicate stage output");
   result=JSON.parse(message);
 }});
 if(result===undefined)throw Error("missing stage output");
 return result;
}

function thirteenScreen(suppliedStageInput,console){
// Run with Node.js. Outputs all remaining explicit target sets in the 13-free branch.
// Node.js standard library only. Exact five non-Euler exceptions with 3 absent.
let data;
{
const X=1e6,P=500;
const is=new Uint8Array(X+1).fill(1,2),rem=new Float64Array(X+1),fs=new Map(),ps=[];
for(let p=2;p*p<=X;p++)if(is[p])for(let j=p*p;j<=X;j+=p)is[j]=0;
for(let p=5;p<=X;p++)if(is[p]){ps.push(p);rem[p]=p*p+p+1;fs.set(p,[]);}
function pow(a,e,p){let b=1;while(e){if(e%2)b=b*a%p;a=a*a%p;e=Math.floor(e/2);}return b;}
for(const p of [3,...ps]){let a=[];if(p===3)a=[1];else if(p%3===1){for(let k=2;k<p;k++){let z=pow(k,(p-1)/3,p);if(z!==1){a=[z,p-1-z];break;}}}for(const z of a)for(let j=z;j<=X;j+=p)if(rem[j]){fs.get(j).push(p);while(rem[j]%p===0)rem[j]/=p;}}
const caps=new Map();
for(const p of ps){let f=fs.get(p);if(rem[p]>1)f.push(rem[p]);let l=p%3===1?p:f.at(-1);if(p>P)(caps.get(l)??caps.set(l,[]).get(l)).push(p);}

const F=(a,b=1)=>[BigInt(a),BigInt(b)], mul=(a,b)=>[a[0]*b[0],a[1]*b[1]], add=(a,b)=>[a[0]*b[1]+a[1]*b[0],a[1]*b[1]], cmp=(a,b)=>a[0]*b[1]-a[1]*b[0], prod=as=>as.reduce(mul,F(1)), hp=p=>F(p,p-1);
function req(v,msg){if(!v)throw Error(msg);}
function froot(p){let l=1000,h=p*1000;while(h-l>1){let k=Math.floor((l+h)/2);if(BigInt(k)**6n<=BigInt(p)*10n**18n)l=k;else h=k;}return F(l,l-1000);}

const weightBound=F(1006,1000), W=new Map([...caps].map(([l,a])=>[l,prod(a.map(hp))]));
for(const weight of W.values())req(cmp(weight,weightBound)<0,"tail group below 1.006");
const small=ps.filter(p=>p<=P), other=small.filter(p=>p%3===2&&p!==5&&p!==11);
const cand=[...new Set([...small.filter(p=>p%3===1),...small.filter(p=>p%3===2).flatMap(p=>fs.get(p))])].sort((a,b)=>a-b);

const specs=[
 ["5 exceptional, 11 regular, 31 absent",[7,13,19],1,[31]],
 ["5 exceptional, 11 regular, 31 present",[7,13,19,31],1,[]],
 ["5 regular, 11 exceptional, 19 absent",[7,13,31],1,[19]],
 ["5 regular, 11 exceptional, 19 present",[7,13,31,19],1,[]],
 ["both regular, 13 absent",[7,19,31],0,[13]],
 ["both regular, 13 present",[7,19,31,13],0,[]]];
const candSet=new Set(cand);
const stats=[];
const topWeights=[...W.entries()].sort((a,b)=>{let d=cmp(b[1],a[1]);return d>0n?1:d<0n?-1:0;}).slice(0,12);
for(const [name,fixedT,e,banned] of specs.slice(4,5)){
 const setup=cand.filter(p=>!fixedT.includes(p)&&!banned.includes(p));const free=6-e-fixedT.length;
 const anonymous=ps.filter(p=>p%3===1&&!candSet.has(p)&&!banned.includes(p)).slice(0,6);
 let n=0,surv=0,best=F(0),worst,remain=[];
 function rec(i,t){
  const T=[...fixedT,...t], S=new Set(T), L=new Set([7,19,31,...T]);
  const reg=other.filter(p=>fs.get(p).every(r=>S.has(r))), opt=other.filter(p=>!reg.includes(p));
  const active=[5,7,11,...reg,...T.filter(p=>p!==7&&p<=P)];
  for(let u=0;u<=free-t.length;u++){
   if(T.length+u<5||(!T.some(p=>p%12===1)&&!u))continue;
   const used=opt.slice(0,free-t.length-u);
   const ub=prod([...active.map(hp),...used.map(hp),...[...L].map(r=>W.get(r)??F(1)),...topWeights.filter(([r])=>!L.has(r)).slice(0,6-L.size).map(([,w])=>w)]);
   const d=add(F(6,X),mul(F(X+2,X*10000),prod([...T,...anonymous.slice(0,u)].map(froot))));
   const threshold=add(ub,mul(F(2),d));
   n++;if(cmp(threshold,F(2))>=0){surv++;remain.push([T,u]);}
   if(cmp(threshold,best)>0){best=threshold;worst=[T,u,active,used,Number(ub[0])/Number(ub[1]),Number(d[0])/Number(d[1])];}
  }
  if(t.length<free)for(let j=i;j<setup.length;j++)rec(j+1,[...t,setup[j]]);
 }
 rec(0,[]);
 stats.push({name,n,surv,threshold:Number(best[0])/Number(best[1]),worst,remain});
} 


const T=[7,19,31,79,307], S=new Set(T), C=new Set(cand);
const remainder=new Map(), common=[];
for(const p of ps)if(p%3===2){
 const missing=fs.get(p).filter(r=>!S.has(r));
 if(missing.length===0)common.push(p);
 else if(missing.length===1&&missing[0]%12===1&&!C.has(missing[0]))(remainder.get(missing[0])??remainder.set(missing[0],[]).get(missing[0])).push(p);
}
const commonU=prod([...T,...common].map(hp));
const d=q=>add(F(1,X),mul(F(X+2,X*10000),prod([...T,q].map(froot))));
const lower=ps.find(p=>p%12===1&&!C.has(p)&&p!==13);
const gap=cmp(add(mul(commonU,hp(lower)),mul(F(2),d(lower))),F(2));
const possible=[...remainder].filter(([q,bases])=>cmp(add(mul(mul(commonU,hp(q)),prod(bases.map(hp))),mul(F(2),d(q))),F(2))>=0);

req(gap<0n,"anonymous Euler with no middle-support prime");
req(possible.length===1&&possible[0][0]===541,"exception to anonymous Euler screen");
const five=stats[0].remain.filter(([t,u])=>t.length+u===5);
const six=stats[0].remain.filter(([t,u])=>t.length+u===6);
req(five.length===24&&five.every(([,u])=>!u),"five-target screen");
const hidden=six.filter(([,u])=>u);
req(six.length===32&&hidden.length===1&&hidden[0][1]===1&&hidden[0][0].join(",")===T.join(","),"six-target screen");
console.log(JSON.stringify({five:five.map(([t])=>t),six:[...six.filter(([,u])=>!u).map(([t])=>t),[...T,541]],middleGroups:remainder.size,common,lower,exception:possible}));
}
}

function thirteenIntervals(suppliedStageInput,console){
// Run with Node.js and pass the target JSON output of the first stage on stdin.
const X=10n**16n, max=X*X+X+1n;
const G=([a,b])=>{let x=a,y=b;while(y)[x,y]=[y,x%y];return [a/x,b/x];};
const F=(a,b=1n)=>[BigInt(a),BigInt(b)], mul=(a,b)=>G([a[0]*b[0],a[1]*b[1]]), cmp=(a,b)=>a[0]*b[1]-a[1]*b[0];
function isqrt(n){if(n<2n)return n;let g=1n<<BigInt(Math.ceil(n.toString(2).length/2));let h=(g+n/g)>>1n;while(h<g){g=h;h=(g+n/g)>>1n;}return g;}
function smooth(T){
 let ans=[];
 function go(i,v){
  if(i===T.length){let d=4n*v-3n,r=isqrt(d),q=(r-1n)/2n;if(r*r===d&&q%6n===5n&&q<=X)ans.push(q);return;}
  for(let w=v;w<=max;w*=BigInt(T[i]))go(i+1,w);
 }
 go(0,1n);return ans;
}
function ints(p,role,optional){
 p=BigInt(p);const v=optional?[[F(1),F(1)]]:[];
 for(let e=role==="Euler"?1:2;;e+=role==="Euler"?4:2){
  if(((e+1)%3===0)!==(role==="reg"))continue;
  let pow=p**BigInt(e),val=F(p*pow-1n,(p-1n)*pow);
  if((p-1n)*pow>10n**10n){v.push([val,F(p,p-1n)]);break;}
  v.push([val,val]);
 }
 return v;
}
if(8n**4n>=7n*5n**4n||(X+2n)*8n**5n>=140n*X*3n**5n)throw Error("fine tail bound");
const hiThreshold=F(2); // 2(1 - 140/10^8)

function req(condition,message){if(!condition)throw Error(message);}


function order(p,r){
 let x=Number(p%BigInt(r)),y=1;
 for(let k=1;k<r;k++){y=y*x%r;if(y===1)return k;}
 throw Error("order");
}
function regular(T,p,optional){
 const es=[...new Set(T.map(r=>order(p,r)).filter(n=>n%6===3))].sort((a,b)=>a-b).map(n=>n-1);
 const choices=optional?[[F(1),F(1)]]:[];
 for(const e of es){
  let pow=p**BigInt(e), h=F(p*pow-1n,(p-1n)*pow);
  if((p-1n)*pow>10n**10n){choices.push([h,F(p,p-1n)]);break;}
  choices.push([h,h]);
 }
 return choices;
}
function prime(p){for(let j=2;j*j<=p;j++)if(p%j===0)return false;return true;}
function boxes(T,pi,R,q,B){
 const count=T.length;
 const fine=count===5?140n:360n;
 if((X+2n)*8n**BigInt(count)>=fine*X*3n**BigInt(count))throw Error("fine tail bound");
 const loThreshold=F(2n*(10n**8n-fine),10n**8n);
 const domains=[
  ...T.map(p=>ints(p,p===pi?"Euler":"ex",false)),
  ...(q?[ints(q,"ex",false)]:[]),
  ...R.filter(p=>BigInt(q)!==p).map(p=>regular(T,p,p!==5n&&p!==11n))
 ];
 if(domains.some(x=>!x.length))return null;
 const suffix=Array(domains.length+1);suffix[domains.length]=q||T.length===6?F(1):F(B+1,B);
 for(let i=domains.length-1;i>=0;i--)suffix[i]=mul(suffix[i+1],domains[i].at(-1)[1]);
 let best=F(0),num=0;
 function go(i,lo,hi){
  if(cmp(lo,hiThreshold)>0||cmp(mul(hi,suffix[i]),loThreshold)<0)return;
  if(i===domains.length){num++;if(cmp(hi,best)>0)best=hi;return;}
  for(const [a,b] of domains[i])go(i+1,mul(lo,a),mul(hi,b));
 }
 go(0,F(1),F(1));if(!num)return null;
 const gap=loThreshold[0]*best[1]-loThreshold[1]*best[0];
 return [q,num,gap>0n?String((loThreshold[0]*best[1])/gap):"unbounded"];
}

const B=500, qs=Array.from({length:B},(_,i)=>i+1).filter(p=>p>11&&p%3===2&&prime(p));
const input=readInput(),out=[];
for(const [k,sets] of Object.entries({five:input.five,six:input.six})){
 let count=0;
 for(const T of sets){
  const R=smooth(T);if(!R.includes(5n)||!R.includes(11n))throw Error("missing required preimage");
  for(const pi of T.filter(p=>p%12===1)){
   count++;
   for(const q of k==="five"?[...qs,0]:[0])if(boxes(T,pi,R,q,B))throw Error("interval survived "+T+" "+pi+" "+q);
  }
 }
 out.push([k,sets.length,count]);
}
console.log(JSON.stringify(out));

function readInput(){return suppliedStageInput;}
}

function anonymousScreen(suppliedStageInput,console){
// Run with Node.js. Check the unlocated target families after forcing 13.
// Run with Node.js. Outputs all remaining explicit target sets in the 13-free branch.
// Node.js standard library only. Exact five non-Euler exceptions with 3 absent.
let data;
{
const X=1e6,P=500;
const is=new Uint8Array(X+1).fill(1,2),rem=new Float64Array(X+1),fs=new Map(),ps=[];
for(let p=2;p*p<=X;p++)if(is[p])for(let j=p*p;j<=X;j+=p)is[j]=0;
for(let p=5;p<=X;p++)if(is[p]){ps.push(p);rem[p]=p*p+p+1;fs.set(p,[]);}
function pow(a,e,p){let b=1;while(e){if(e%2)b=b*a%p;a=a*a%p;e=Math.floor(e/2);}return b;}
for(const p of [3,...ps]){let a=[];if(p===3)a=[1];else if(p%3===1){for(let k=2;k<p;k++){let z=pow(k,(p-1)/3,p);if(z!==1){a=[z,p-1-z];break;}}}for(const z of a)for(let j=z;j<=X;j+=p)if(rem[j]){fs.get(j).push(p);while(rem[j]%p===0)rem[j]/=p;}}
const caps=new Map();
for(const p of ps){let f=fs.get(p);if(rem[p]>1)f.push(rem[p]);let l=p%3===1?p:f.at(-1);if(p>P)(caps.get(l)??caps.set(l,[]).get(l)).push(p);}

const F=(a,b=1)=>[BigInt(a),BigInt(b)], mul=(a,b)=>[a[0]*b[0],a[1]*b[1]], add=(a,b)=>[a[0]*b[1]+a[1]*b[0],a[1]*b[1]], cmp=(a,b)=>a[0]*b[1]-a[1]*b[0], prod=as=>as.reduce(mul,F(1)), hp=p=>F(p,p-1);
function req(v,msg){if(!v)throw Error(msg);}
function froot(p){let l=1000,h=p*1000;while(h-l>1){let k=Math.floor((l+h)/2);if(BigInt(k)**6n<=BigInt(p)*10n**18n)l=k;else h=k;}return F(l,l-1000);}

const weightBound=F(1006,1000), W=new Map([...caps].map(([l,a])=>[l,prod(a.map(hp))]));
for(const weight of W.values())req(cmp(weight,weightBound)<0,"tail group below 1.006");
const small=ps.filter(p=>p<=P), other=small.filter(p=>p%3===2&&p!==5&&p!==11);
const cand=[...new Set([...small.filter(p=>p%3===1),...small.filter(p=>p%3===2).flatMap(p=>fs.get(p))])].sort((a,b)=>a-b);

const specs=[
 ["5 exceptional, 11 regular, 31 absent",[7,13,19],1,[31]],
 ["5 exceptional, 11 regular, 31 present",[7,13,19,31],1,[]],
 ["5 regular, 11 exceptional, 19 absent",[7,13,31],1,[19]],
 ["5 regular, 11 exceptional, 19 present",[7,13,31,19],1,[]],
 ["both regular, 13 absent",[7,19,31],0,[13]],
 ["both regular, 13 present",[7,19,31,13],0,[]]];
const candSet=new Set(cand);
const stats=[];
const topWeights=[...W.entries()].sort((a,b)=>{let d=cmp(b[1],a[1]);return d>0n?1:d<0n?-1:0;}).slice(0,12);
for(const [name,fixedT,e,banned] of specs.slice(5)){
 const setup=cand.filter(p=>!fixedT.includes(p)&&!banned.includes(p));const free=6-e-fixedT.length;
 const anonymous=ps.filter(p=>p%3===1&&!candSet.has(p)&&!banned.includes(p)).slice(0,6);
 let n=0,surv=0,best=F(0),worst,remain=[];
 function rec(i,t){
  const T=[...fixedT,...t], S=new Set(T), L=new Set([7,19,31,...T]);
  const reg=other.filter(p=>fs.get(p).every(r=>S.has(r))), opt=other.filter(p=>!reg.includes(p));
  const active=[5,7,11,...reg,...T.filter(p=>p!==7&&p<=P)];
  for(let u=0;u<=free-t.length;u++){
   if(T.length+u<5||(!T.some(p=>p%12===1)&&!u))continue;
   const used=opt.slice(0,free-t.length-u);
   const ub=prod([...active.map(hp),...used.map(hp),...[...L].map(r=>W.get(r)??F(1)),...topWeights.filter(([r])=>!L.has(r)).slice(0,6-L.size).map(([,w])=>w)]);
   const d=add(F(6,X),mul(F(X+2,X*10000),prod([...T,...anonymous.slice(0,u)].map(froot))));
   const threshold=add(ub,mul(F(2),d));
   n++;if(cmp(threshold,F(2))>=0){surv++;remain.push([T,u]);}
   if(cmp(threshold,best)>0){best=threshold;worst=[T,u,active,used,Number(ub[0])/Number(ub[1]),Number(d[0])/Number(d[1])];}
  }
  if(t.length<free)for(let j=i;j<setup.length;j++)rec(j+1,[...t,setup[j]]);
 }
 rec(0,[]);
 stats.push({name,n,surv,threshold:Number(best[0])/Number(best[1]),worst,remain});
} 



const hidden=stats[0].remain.filter(([T,u])=>u);
req(hidden.length===6&&hidden.every(([,u])=>u===1),"unlocated target families");
const out=[];
for(const [T] of hidden){
 const C=new Set(cand),S=new Set(T),remainder=new Map(),common=[];
 for(const p of ps)if(p%3===2){
  const missing=fs.get(p).filter(r=>!S.has(r));
  if(missing.length===0)common.push(p);
  else if(missing.length===1&&!C.has(missing[0]))(remainder.get(missing[0])??remainder.set(missing[0],[]).get(missing[0])).push(p);
 }
 const commonU=prod([...T,...common].map(hp));
 
 const lower=ps.find(p=>p%3===1&&!C.has(p));
 const scale=10n**15n, Y2=10n**18n;
 let acc=0n, terms=0;
 function isqrt(n){if(n<2n)return n;let x=1n<<BigInt(Math.ceil(n.toString(2).length/2));let y=(x+n/x)>>1n;while(y<x){x=y;y=(x+n/x)>>1n;}return x;}
 const bases=[...T,lower].map(BigInt);
 function walk(i,m){
  if(i===bases.length){let root=isqrt(m), div=root>BigInt(X)?root:BigInt(X);acc+=scale/div+1n;terms++;return;}
  for(let n=m;n<=Y2;n*=bases[i])walk(i+1,n);
 }
 walk(0,1n);
 const tailKnown=add(F(6-T.length,X),mul(F(X+2,X),add(F(acc,scale),mul(F(1,1000000),prod([...T,lower].map(froot))))));
 const d=q=>tailKnown;
 const u=(q,p)=>add(mul(mul(commonU,hp(q)),prod(p.map(hp))),mul(F(2),d(q)));
 const possible=[...remainder].filter(([q,bases])=>cmp(u(q,bases),F(2))>=0);
 let check=u(lower,[]);
 let upper=mul(commonU,hp(lower));
 for(const [p,ps] of remainder){const u=mul(mul(commonU,hp(p)),prod(ps.map(hp)));if(cmp(u,upper)>0n)upper=u;}
 const limit=add(F(2),mul(F(-2),tailKnown));const diff=cmp(limit,upper);
 const bMax=diff>0n?(limit[0]*upper[1])/diff:0n;
 if(T.length===4)req(bMax===24n,"bound on the remaining exceptional base");
 if(T.length===5&&(T.includes(37)||T.includes(43)))req(cmp(check,F(2))<0n&&possible.length===0,"anonymous six-target exclusion");
 out.push({T,bMax:T.length===4?String(bMax):null,noExtraExceptionUpper:upper.map(String),lower,common,empty:Number(check[0])/Number(check[1]),sign:String(cmp(check,F(2))),count:remainder.size,remaining:possible.length,possible:possible.slice(0,12),terms,D:Number(tailKnown[0])/Number(tailKnown[1])});
}
console.log(JSON.stringify(out));}
}

function excludeTwentyThree(suppliedStageInput,console){
// Run with Node. Both sections use exact integer comparisons.
let candidates;
{
const console={log:s=>{candidates=JSON.parse(s).rs;}};
// Run with Node.js. Check the unlocated target families after forcing 13.
// Run with Node.js. Outputs all remaining explicit target sets in the 13-free branch.
// Node.js standard library only. Exact five non-Euler exceptions with 3 absent.
let data;
{
const X=1e6,P=500;
const is=new Uint8Array(X+1).fill(1,2),rem=new Float64Array(X+1),fs=new Map(),ps=[];
for(let p=2;p*p<=X;p++)if(is[p])for(let j=p*p;j<=X;j+=p)is[j]=0;
for(let p=5;p<=X;p++)if(is[p]){ps.push(p);rem[p]=p*p+p+1;fs.set(p,[]);}
function pow(a,e,p){let b=1;while(e){if(e%2)b=b*a%p;a=a*a%p;e=Math.floor(e/2);}return b;}
for(const p of [3,...ps]){let a=[];if(p===3)a=[1];else if(p%3===1){for(let k=2;k<p;k++){let z=pow(k,(p-1)/3,p);if(z!==1){a=[z,p-1-z];break;}}}for(const z of a)for(let j=z;j<=X;j+=p)if(rem[j]){fs.get(j).push(p);while(rem[j]%p===0)rem[j]/=p;}}
const caps=new Map();
for(const p of ps){let f=fs.get(p);if(rem[p]>1)f.push(rem[p]);let l=p%3===1?p:f.at(-1);if(p>P)(caps.get(l)??caps.set(l,[]).get(l)).push(p);}

const F=(a,b=1)=>[BigInt(a),BigInt(b)], mul=(a,b)=>[a[0]*b[0],a[1]*b[1]], add=(a,b)=>[a[0]*b[1]+a[1]*b[0],a[1]*b[1]], cmp=(a,b)=>a[0]*b[1]-a[1]*b[0], prod=as=>as.reduce(mul,F(1)), hp=p=>F(p,p-1);
function req(v,msg){if(!v)throw Error(msg);}
function froot(p){let l=1000,h=p*1000;while(h-l>1){let k=Math.floor((l+h)/2);if(BigInt(k)**6n<=BigInt(p)*10n**18n)l=k;else h=k;}return F(l,l-1000);}

const weightBound=F(1006,1000), W=new Map([...caps].map(([l,a])=>[l,prod(a.map(hp))]));
for(const weight of W.values())req(cmp(weight,weightBound)<0,"tail group below 1.006");
const small=ps.filter(p=>p<=P), other=small.filter(p=>p%3===2&&p!==5&&p!==11);
const cand=[...new Set([...small.filter(p=>p%3===1),...small.filter(p=>p%3===2).flatMap(p=>fs.get(p))])].sort((a,b)=>a-b);

const specs=[
 ["5 exceptional, 11 regular, 31 absent",[7,13,19],1,[31]],
 ["5 exceptional, 11 regular, 31 present",[7,13,19,31],1,[]],
 ["5 regular, 11 exceptional, 19 absent",[7,13,31],1,[19]],
 ["5 regular, 11 exceptional, 19 present",[7,13,31,19],1,[]],
 ["both regular, 13 absent",[7,19,31],0,[13]],
 ["both regular, 13 present",[7,19,31,13],0,[]]];
const candSet=new Set(cand);
const stats=[];
const topWeights=[...W.entries()].sort((a,b)=>{let d=cmp(b[1],a[1]);return d>0n?1:d<0n?-1:0;}).slice(0,12);
for(const [name,fixedT,e,banned] of specs.slice(5)){
 const setup=cand.filter(p=>!fixedT.includes(p)&&!banned.includes(p));const free=6-e-fixedT.length;
 const anonymous=ps.filter(p=>p%3===1&&!candSet.has(p)&&!banned.includes(p)).slice(0,6);
 let n=0,surv=0,best=F(0),worst,remain=[];
 function rec(i,t){
  const T=[...fixedT,...t], S=new Set(T), L=new Set([7,19,31,...T]);
  const reg=other.filter(p=>fs.get(p).every(r=>S.has(r))), opt=other.filter(p=>!reg.includes(p));
  const active=[5,7,11,...reg,...T.filter(p=>p!==7&&p<=P)];
  for(let u=0;u<=free-t.length;u++){
   if(T.length+u<5||(!T.some(p=>p%12===1)&&!u))continue;
   const used=opt.slice(0,free-t.length-u);
   const ub=prod([...active.map(hp),...used.map(hp),...[...L].map(r=>W.get(r)??F(1)),...topWeights.filter(([r])=>!L.has(r)).slice(0,6-L.size).map(([,w])=>w)]);
   const d=add(F(6,X),mul(F(X+2,X*10000),prod([...T,...anonymous.slice(0,u)].map(froot))));
   const threshold=add(ub,mul(F(2),d));
   n++;if(cmp(threshold,F(2))>=0){surv++;remain.push([T,u]);}
   if(cmp(threshold,best)>0){best=threshold;worst=[T,u,active,used,Number(ub[0])/Number(ub[1]),Number(d[0])/Number(d[1])];}
  }
  if(t.length<free)for(let j=i;j<setup.length;j++)rec(j+1,[...t,setup[j]]);
 }
 rec(0,[]);
 stats.push({name,n,surv,threshold:Number(best[0])/Number(best[1]),worst,remain});
} 



const hidden=stats[0].remain.filter(([T,u])=>u);
req(hidden.length===6&&hidden.every(([,u])=>u===1),"unlocated target families");
const out=[];
for(const [T] of hidden){
 const C=new Set(cand),S=new Set(T),remainder=new Map(),common=[];
 for(const p of ps)if(p%3===2){
  const missing=fs.get(p).filter(r=>!S.has(r));
  if(missing.length===0)common.push(p);
  else if(missing.length===1&&!C.has(missing[0]))(remainder.get(missing[0])??remainder.set(missing[0],[]).get(missing[0])).push(p);
 }
 const commonU=prod([...T,...common].map(hp));
 
 const lower=ps.find(p=>p%3===1&&!C.has(p));
 const scale=10n**15n, Y2=10n**18n;
 let acc=0n, terms=0;
 function isqrt(n){if(n<2n)return n;let x=1n<<BigInt(Math.ceil(n.toString(2).length/2));let y=(x+n/x)>>1n;while(y<x){x=y;y=(x+n/x)>>1n;}return x;}
 const bases=[...T,lower].map(BigInt);
 function walk(i,m){
  if(i===bases.length){let root=isqrt(m), div=root>BigInt(X)?root:BigInt(X);acc+=scale/div+1n;terms++;return;}
  for(let n=m;n<=Y2;n*=bases[i])walk(i+1,n);
 }
 walk(0,1n);
 const tailKnown=add(F(6-T.length,X),mul(F(X+2,X),add(F(acc,scale),mul(F(1,1000000),prod([...T,lower].map(froot))))));
 const d=q=>tailKnown;
 const u=(q,p)=>add(mul(mul(commonU,hp(q)),prod(p.map(hp))),mul(F(2),d(q)));
 const possible=[...remainder].filter(([q,bases])=>cmp(u(q,bases),F(2))>=0);
 let check=u(lower,[]);
 if(T.length===4){
 const lim=add(F(2),mul(F(-2),tailKnown)), b=23;
 const c=mul(commonU,hp(b)),g=cmp(lim,c);
 req(g>0n && lim[0]*c[1]/g===5414n,"uniform empty-group bound for b=23");
 const test=(q,bases)=>cmp(add(mul(mul(c,hp(q)),prod(bases.map(hp))),mul(F(2),tailKnown)),F(2))>=0;
 const rs=[...new Set([...ps.filter(q=>q<=5414&&q%3===1&&!C.has(q)&&test(q,remainder.get(q)??[])),...[...remainder].filter(([q,bases])=>test(q,bases)).map(([q])=>q)])].sort((a,b)=>a-b);
 out.push({rs});
 }
 let upper=mul(commonU,hp(lower));
 for(const [p,ps] of remainder){const u=mul(mul(commonU,hp(p)),prod(ps.map(hp)));if(cmp(u,upper)>0n)upper=u;}
 const limit=add(F(2),mul(F(-2),tailKnown));const diff=cmp(limit,upper);
 const bMax=diff>0n?(limit[0]*upper[1])/diff:0n;
 if(T.length===4)req(bMax===24n,"bound on the remaining exceptional base");
 if(T.length===5&&(T.includes(37)||T.includes(43)))req(cmp(check,F(2))<0n&&possible.length===0,"anonymous six-target exclusion");
 out.push({T,bMax:T.length===4?String(bMax):null,noExtraExceptionUpper:upper.map(String),lower,common,empty:Number(check[0])/Number(check[1]),sign:String(cmp(check,F(2))),count:remainder.size,remaining:possible.length,possible:possible.slice(0,12),terms,D:Number(tailKnown[0])/Number(tailKnown[1])});
}
console.log(JSON.stringify(out[0]));}
}
{
const rs=candidates;
function put(o){if(o.n!==653||o.out.length!==0)throw Error("unexpected survivor/count");console.log(JSON.stringify({candidateTargets:rs.length,eulerChoices:o.n,excluded:o.out.length===0}));}
if(rs.length!==440)throw Error("candidate-target count");
// Run with Node.js and pass the target JSON output of the first stage on stdin.
const X=10n**16n, max=X*X+X+1n;
const G=([a,b])=>{let x=a,y=b;while(y)[x,y]=[y,x%y];return [a/x,b/x];};
const F=(a,b=1n)=>[BigInt(a),BigInt(b)], mul=(a,b)=>G([a[0]*b[0],a[1]*b[1]]), cmp=(a,b)=>a[0]*b[1]-a[1]*b[0];
function isqrt(n){if(n<2n)return n;let g=1n<<BigInt(Math.ceil(n.toString(2).length/2));let h=(g+n/g)>>1n;while(h<g){g=h;h=(g+n/g)>>1n;}return g;}
const ok=Array.from({length:576},(_,i)=>i).reduce((s,i)=>(s.add(i*i%576),s),new Set());
function smooth(T){
 let ans=[];
 function go(i,v){
  if(i===T.length){let d=4n*v-3n;if(!ok.has(Number(d%576n)))return;let r=isqrt(d),q=(r-1n)/2n;if(r*r===d&&q%6n===5n&&q<=X)ans.push(q);return;}
  for(let w=v;w<=max;w*=BigInt(T[i]))go(i+1,w);
 }
 go(0,1n);return ans;
}
function ints(p,role,optional){
 p=BigInt(p);const v=optional?[[F(1),F(1)]]:[];
 for(let e=role==="Euler"?1:2;;e+=role==="Euler"?4:2){
  if(((e+1)%3===0)!==(role==="reg"))continue;
  let pow=p**BigInt(e),val=F(p*pow-1n,(p-1n)*pow);
  if((p-1n)*pow>10n**10n){v.push([val,F(p,p-1n)]);break;}
  v.push([val,val]);
 }
 return v;
}
if(8n**4n>=7n*5n**4n||(X+2n)*8n**5n>=140n*X*3n**5n)throw Error("fine tail bound");
const hiThreshold=F(2); // 2(1 - 140/10^8)

function req(condition,message){if(!condition)throw Error(message);}


function order(p,R){
 const r=BigInt(R), x=p%r; if(!x)throw Error("nonunit");
 const pow=(e)=>{let a=x,b=1n;for(let n=BigInt(e);n;n>>=1n,a=a*a%r)if(n&1n)b=b*a%r;return b};
 let ord=R-1,k=ord;
 for(let d=2;d*d<=k;d++)if(k%d===0){while(k%d===0)k/=d;while(ord%d===0&&pow(ord/d)===1n)ord/=d;}
 if(k>1)while(ord%k===0&&pow(ord/k)===1n)ord/=k;
 return ord;
}
function regular(T,p,optional){
 const es=[...new Set(T.map(r=>order(p,r)).filter(n=>n%6===3))].sort((a,b)=>a-b).map(n=>n-1);
 const choices=optional?[[F(1),F(1)]]:[];
 for(const e of es){
  let pow=p**BigInt(e), h=F(p*pow-1n,(p-1n)*pow);
  if((p-1n)*pow>10n**10n){choices.push([h,F(p,p-1n)]);break;}
  choices.push([h,h]);
 }
 return choices;
}
function prime(p){for(let j=2;j*j<=p;j++)if(p%j===0)return false;return true;}
function boxes(T,pi,R,q,B){
 const count=T.length;
 const fine=count===5?140n:360n;
 if((X+2n)*8n**BigInt(count)>=fine*X*3n**BigInt(count))throw Error("fine tail bound");
 const loThreshold=F(2n*(10n**8n-fine),10n**8n);
 const domains=[
  ...T.map(p=>ints(p,p===pi?"Euler":"ex",false)),
  ...(q?[ints(q,"ex",false)]:[]),
  ...R.filter(p=>BigInt(q)!==p).map(p=>regular(T,p,p!==5n&&p!==11n))
 ];
 if(domains.some(x=>!x.length))return null;
 const suffix=Array(domains.length+1);suffix[domains.length]=q||T.length===6?F(1):F(B+1,B);
 for(let i=domains.length-1;i>=0;i--)suffix[i]=mul(suffix[i+1],domains[i].at(-1)[1]);
 let best=F(0),num=0;
 function go(i,lo,hi){if(num)return;
  if(cmp(lo,hiThreshold)>0||cmp(mul(hi,suffix[i]),loThreshold)<0)return;
  if(i===domains.length){num++;if(cmp(hi,best)>0)best=hi;return;}
  for(const [a,b] of domains[i])go(i+1,mul(lo,a),mul(hi,b));
 }
 go(0,F(1),F(1));if(!num)return null;
 const gap=loThreshold[0]*best[1]-loThreshold[1]*best[0];
 return [q,num,gap>0n?String((loThreshold[0]*best[1])/gap):"unbounded"];
}

const B=500, qs=Array.from({length:B},(_,i)=>i+1).filter(p=>p>11&&p%3===2&&prime(p));

let n=0,out=[];
for(const r of rs){
 const T=[7,13,19,31,r],R=smooth(T).filter(p=>T.every(r=>p%BigInt(r)!==0n));
 if(!R.includes(5n)||!R.includes(11n))throw Error("mandatory base");
 for(const pi of T.filter(p=>p%12===1)){
   n++;
   const v=boxes(T,pi,R,23,500);
   if(v)out.push({r,pi,R:R.map(String)});
 }
}
put({n,out});

}
}

function seventeenScreen(suppliedStageInput,console){
// Run with Node.js. Outputs exact anonymous-five b=17 candidates for the interval checker.
{
function save(stats){if(stats.length!==2||stats.some(([_,best,value,cross])=>best.join(",")!=="2453195303377,1226928384000"||cross))throw Error("core abundance maximum");}
const X=10n**16n, max=X*X+X+1n;
const G=([a,b])=>{let x=a,y=b;while(y)[x,y]=[y,x%y];return [a/x,b/x];};
const F=(a,b=1n)=>[BigInt(a),BigInt(b)], mul=(a,b)=>G([a[0]*b[0],a[1]*b[1]]), cmp=(a,b)=>a[0]*b[1]-a[1]*b[0];
function isqrt(n){if(n<2n)return n;let g=1n<<BigInt(Math.ceil(n.toString(2).length/2));let h=(g+n/g)>>1n;while(h<g){g=h;h=(g+n/g)>>1n;}return g;}
const ok=Array.from({length:576},(_,i)=>i).reduce((s,i)=>(s.add(i*i%576),s),new Set());
function smooth(T){
 let ans=[];
 function go(i,v){
  if(i===T.length){let d=4n*v-3n;if(!ok.has(Number(d%576n)))return;let r=isqrt(d),q=(r-1n)/2n;if(r*r===d&&q%6n===5n&&q<=X)ans.push(q);return;}
  for(let w=v;w<=max;w*=BigInt(T[i]))go(i+1,w);
 }
 go(0,1n);return ans;
}
function ints(p,role,optional){
 p=BigInt(p);const v=optional?[[F(1),F(1),0,false]]:[];
 for(let e=role==="Euler"?1:2;;e+=role==="Euler"?4:2){
  if(((e+1)%3===0)!==(role==="reg"))continue;
  let pow=p**BigInt(e),val=F(p*pow-1n,(p-1n)*pow);
  if((p-1n)*pow>10n**10n){v.push([val,F(p,p-1n),e,true]);break;}
  v.push([val,val,e,false]);
 }
 return v;
}
if(8n**4n>=7n*5n**4n||(X+2n)*8n**5n>=140n*X*3n**5n)throw Error("fine tail bound");

const stats=[];
const R=smooth([7,13,19,31]).map(Number).sort((a,b)=>a-b);
if(R.join(',')!=='5,11,191,653')throw Error('core smooth preimages');
for(const externalEuler of [false,true]){
 const rows=[5,7,11,13,17,19,31,191,653].map(p=>[p,ints(p,p===13&&!externalEuler?"Euler":[5,11,191,653].includes(p)?"reg":"ex",p>31)]);
 const n=rows.length, suffix=Array(n+1);suffix[n]=F(1);
 for(let i=n-1;i>=0;i--)suffix[i]=mul(suffix[i+1],rows[i][1].at(-1)[1]);
 let best=F(0),cross=[];
 function go(i,l,h,es){if(cmp(l,F(2))>0||cmp(mul(h,suffix[i]),F(199,100))<0)return;
  if(i===n){
   if(es[3][0]===1&&([7,8].some(j=>es[j][0])))return;
   if(cmp(h,best)>0)best=h;
   if(cmp(h,F(2))>=0)cross.push(es);
   return;
  }
  for(const [a,b,e,t] of rows[i][1])go(i+1,mul(l,a),mul(h,b),[...es,[e,t]]);
 }
 go(0,F(1),F(1),[]);stats.push([externalEuler,best.map(String),Number(best[0])/Number(best[1]),cross.length,cross.slice(0,3)]);
}
save(stats);
}
const lim=1e6,P=500;
const is=new Uint8Array(lim+1).fill(1,2),rem=new Float64Array(lim+1),fs=new Map(),ps=[];
for(let p=2;p*p<=lim;p++)if(is[p])for(let j=p*p;j<=lim;j+=p)is[j]=0;
for(let p=5;p<=lim;p++)if(is[p]){ps.push(p);rem[p]=p*p+p+1;fs.set(p,[]);}
function pow(a,e,p){let b=1;while(e){if(e%2)b=b*a%p;a=a*a%p;e=Math.floor(e/2);}return b;}
for(const p of [3,...ps]){let a=[];if(p===3)a=[1];else if(p%3===1){for(let k=2;k<p;k++){let z=pow(k,(p-1)/3,p);if(z!==1){a=[z,p-1-z];break;}}}for(const z of a)for(let j=z;j<=lim;j+=p)if(rem[j]){fs.get(j).push(p);while(rem[j]%p===0)rem[j]/=p;}}

const X=10n**6n,Y2=10n**18n,M=10n**15n,S=[7n,13n,19n,31n];
const F=(a,b=1n)=>[BigInt(a),BigInt(b)],mul=(a,b)=>[a[0]*b[0],a[1]*b[1]],add=(a,b)=>[a[0]*b[1]+a[1]*b[0],a[1]*b[1]],cmp=(a,b)=>a[0]*b[1]-a[1]*b[0],prod=as=>as.reduce(mul,F(1));
function isqrt(n){if(n<2n)return n;let a=1n<<BigInt(Math.ceil(n.toString(2).length/2)),b=(a+n/a)>>1n;while(b<a){a=b;b=(a+n/a)>>1n;}return a;}
function sixth(n){let a=1000n,b=1000n*n;while(b-a>1n){let c=(a+b)>>1n;if(c**6n<=n*10n**18n)a=c;else b=c;}return a;}
const baseProd=prod(S.map(p=>{let z=sixth(p);return F(z,z-1000n);}));
const H=F(2453195303377n,1226928384000n);
function tail(q){
 let sum=0n,n=0;
 function walk(i,m){if(i===S.length){let d=isqrt(m);if(d<X)d=X;sum+=M/d+1n;n++;return;}
 for(let w=m;w<=Y2;w*=S[i])walk(i+1,w);}
 for(let w=q;w<=Y2;w*=q)walk(0,w);
 const d=mul(F(X+2n,X),add(F(sum,M),mul(F(1,1000000),mul(baseProd,F(1000n,sixth(q)-1000n)))));
 return [add(d,F(140,10n**8n)),n];
}
function pass(q){let [d]=tail(q);return cmp(add(mul(H,F(q,q-1n)),mul(F(2),d)),F(2))<0;}

let C=new Set([...ps.filter(p=>p<=500&&p%3===1),...ps.filter(p=>p<=500&&p%3===2).flatMap(p=>{let f=fs.get(p);return rem[p]>1?[...f,rem[p]]:f;})]);
function eligible(r){return !C.has(r)&&(r%51===1||r%204===169);}
const coreSet=new Set(S.map(Number)), groups=new Map(), normalized=[];
for(const p of ps)if(p%3===2){
 let f=fs.get(p);if(rem[p]>1)f.push(rem[p]);const missing=f.filter(q=>!coreSet.has(q));
 if(!missing.length)normalized.push(p);
 if(missing.length===1&&eligible(missing[0]))(groups.get(missing[0])??groups.set(missing[0],[]).get(missing[0])).push(p);
}
if(normalized.join(",")!=="5,11,191,653")throw Error("core regular profile");
if(!pass(163932n)||pass(163931n))throw Error("empty group boundary");
let pre=[...new Set([...ps.filter(p=>p<163932&&eligible(p)),...groups.keys()])].sort((a,b)=>a-b);
const keep=[];
for(const p of pre){
 const q=BigInt(p),g=groups.get(p)??[];
 const U=mul(mul(H,F(q,q-1n)),prod(g.map(n=>F(n,n-1))));
 const [d]=tail(q);
 if(cmp(add(U,mul(F(2),d)),F(2))>=0n)keep.push([p,g.length]);
}
if(keep.length!==699||groups.size!==1069||pre.length!==1745)throw Error("screen counts");
console.log(JSON.stringify({num:keep.length,groups:groups.size,pre:pre.length,keep}));
}

function seventeenIntervals(suppliedStageInput,console){
// Run with Node.js, taking the exact candidate-screen JSON from stdin.
const X=10n**16n, max=X*X+X+1n;
const G=([a,b])=>{let x=a,y=b;while(y)[x,y]=[y,x%y];return [a/x,b/x];};
const F=(a,b=1n)=>[BigInt(a),BigInt(b)], mul=(a,b)=>G([a[0]*b[0],a[1]*b[1]]), cmp=(a,b)=>a[0]*b[1]-a[1]*b[0];
function isqrt(n){if(n<2n)return n;let g=1n<<BigInt(Math.ceil(n.toString(2).length/2));let h=(g+n/g)>>1n;while(h<g){g=h;h=(g+n/g)>>1n;}return g;}
const ok=Array.from({length:576},(_,i)=>i).reduce((s,i)=>(s.add(i*i%576),s),new Set());
function smooth(T){
 let ans=[];
 function go(i,v){
  if(i===T.length){let d=4n*v-3n;if(!ok.has(Number(d%576n)))return;let r=isqrt(d),q=(r-1n)/2n;if(r*r===d&&q%6n===5n&&q<=X)ans.push(q);return;}
  for(let w=v;w<=max;w*=BigInt(T[i]))go(i+1,w);
 }
 go(0,1n);return ans;
}
function ints(p,role,optional){
 p=BigInt(p);const v=optional?[[F(1),F(1)]]:[];
 for(let e=role==="Euler"?1:2;;e+=role==="Euler"?4:2){
  if(((e+1)%3===0)!==(role==="reg"))continue;
  let pow=p**BigInt(e),val=F(p*pow-1n,(p-1n)*pow);
  if((p-1n)*pow>10n**10n){v.push([val,F(p,p-1n)]);break;}
  v.push([val,val]);
 }
 return v;
}
if(8n**4n>=7n*5n**4n||(X+2n)*8n**5n>=140n*X*3n**5n)throw Error("fine tail bound");
const hiThreshold=F(2); // 2(1 - 140/10^8)

function req(condition,message){if(!condition)throw Error(message);}


function order(p,R){
 const r=BigInt(R), x=p%r; if(!x)throw Error("nonunit");
 const pow=(e)=>{let a=x,b=1n;for(let n=BigInt(e);n;n>>=1n,a=a*a%r)if(n&1n)b=b*a%r;return b};
 let ord=R-1,k=ord;
 for(let d=2;d*d<=k;d++)if(k%d===0){while(k%d===0)k/=d;while(ord%d===0&&pow(ord/d)===1n)ord/=d;}
 if(k>1)while(ord%k===0&&pow(ord/k)===1n)ord/=k;
 return ord;
}

function regular(T,p,optional){
 const es=[...new Set(T.map(r=>order(p,r)).filter(n=>n%6===3))].sort((a,b)=>a-b).map(n=>n-1);
 const choices=optional?[[F(1),F(1)]]:[];
 let cap=0,maxpow=1n;
 while((p-1n)*maxpow<=10n**10n){cap++;maxpow*=p;}
 for(const e of es){
  let power=p**BigInt(Math.min(e,cap)), h=F(p*power-1n,(p-1n)*power);
  if(e>=cap){choices.push([h,F(p,p-1n)]);break;}
  choices.push([h,h]);
 }
 return choices;
}
function prime(p){for(let j=2;j*j<=p;j++)if(p%j===0)return false;return true;}
function boxes(T,pi,R,q,B){
 const count=T.length;
 const fine=count===5?140n:360n;
 if((X+2n)*8n**BigInt(count)>=fine*X*3n**BigInt(count))throw Error("fine tail bound");
 const loThreshold=F(2n*(10n**8n-fine),10n**8n);
 const domains=[
  ...T.map(p=>ints(p,p===pi?"Euler":"ex",false)),
  ...(q?[ints(q,"ex",false)]:[]),
  ...R.filter(p=>BigInt(q)!==p).map(p=>regular(T,p,p!==5n&&p!==11n))
 ];
 if(domains.some(x=>!x.length))return null;
 const lower=Array(domains.length+1);lower[domains.length]=F(1);
 for(let i=domains.length-1;i>=0;i--)lower[i]=mul(lower[i+1],domains[i][0][0]);
 const suffix=Array(domains.length+1);suffix[domains.length]=q||T.length===6?F(1):F(B+1,B);
 for(let i=domains.length-1;i>=0;i--)suffix[i]=mul(suffix[i+1],domains[i].at(-1)[1]);
 let best=F(0),num=0;
 function go(i,lo,hi){if(num)return;
  if(cmp(mul(lo,lower[i]),hiThreshold)>0||cmp(mul(hi,suffix[i]),loThreshold)<0)return;
  if(i===domains.length){num++;if(cmp(hi,best)>0)best=hi;return;}
  for(const [a,b] of domains[i])go(i+1,mul(lo,a),mul(hi,b));
 }
 go(0,F(1),F(1));if(!num)return null;
 const gap=loThreshold[0]*best[1]-loThreshold[1]*best[0];
 return [q,num,gap>0n?String((loThreshold[0]*best[1])/gap):"unbounded"];
}


let n=0,out=[];
for(const [r] of readInput().keep){
 const T=[7,13,19,31,r],R=smooth(T).filter(p=>T.every(r=>p%BigInt(r)!==0n));
 if(!R.includes(5n)||!R.includes(11n))throw Error("mandatory base");
 const choices=r%17===16?[r]:T.filter(p=>p%12===1);
 for(const pi of choices){
   if(pi%12!==1)throw Error("Euler residue");
   n++;
   if(boxes(T,pi,R,17,500))out.push({r,pi,R:R.map(String)});
 }
}

if(n!==927||out.map(x=>x.r+","+x.pi).join(";")!=="3673,3673;13479913,13479913")throw Error("unexpected interval result");
function check(ok,why){if(!ok)throw Error(why);}
function isPrime(p){for(let d=2;d*d<=p;d++)if(p%d===0)return false;return true;}
check(3674===2*11*167&&isPrime(167)&&167%3===2,"forced base 167");
check(167*167+167+1===28057&&isPrime(28057)&&28057%3===1,"forced missing target");
check(13479914%37===0&&isPrime(37)&&37%3===1,"forced target 37");
console.log(JSON.stringify({eulerChoices:n,out,excluded:true}));
function readInput(){return suppliedStageInput;}
}

function sixThresholds(suppliedStageInput,console){
// Run with Node.js. Bounds the missing prime in the three six-target empty groups.
function check(ok,message){if(!ok)throw Error(message);}
function prime(p){for(let d=2;d*d<=p;d++)if(p%d===0)return false;return true;}
check(2801===1+7+7**2+7**3+7**4&&prime(2801)&&2801%3===2,"forced regular base");
check(2801**2+2801+1===37*43*4933&&[37,43,4933].every(prime),"three missing targets");
function coreMaxima(){
const X=10n**16n, max=X*X+X+1n;
const G=([a,b])=>{let x=a,y=b;while(y)[x,y]=[y,x%y];return [a/x,b/x];};
const F=(a,b=1n)=>[BigInt(a),BigInt(b)], mul=(a,b)=>G([a[0]*b[0],a[1]*b[1]]), cmp=(a,b)=>a[0]*b[1]-a[1]*b[0];
function isqrt(n){if(n<2n)return n;let g=1n<<BigInt(Math.ceil(n.toString(2).length/2));let h=(g+n/g)>>1n;while(h<g){g=h;h=(g+n/g)>>1n;}return g;}
const ok=Array.from({length:576},(_,i)=>i).reduce((s,i)=>(s.add(i*i%576),s),new Set());
function smooth(T){
 let ans=[];
 function go(i,v){
  if(i===T.length){let d=4n*v-3n;if(!ok.has(Number(d%576n)))return;let r=isqrt(d),q=(r-1n)/2n;if(r*r===d&&q%6n===5n&&q<=X)ans.push(q);return;}
  for(let w=v;w<=max;w*=BigInt(T[i]))go(i+1,w);
 }
 go(0,1n);return ans;
}
function ints(p,role,optional){
 p=BigInt(p);const v=optional?[[F(1),F(1),0,false]]:[];
 for(let e=role==="Euler"?1:2;;e+=role==="Euler"?4:2){
  if(((e+1)%3===0)!==(role==="reg"))continue;
  let pow=p**BigInt(e),val=F(p*pow-1n,(p-1n)*pow);
  if((p-1n)*pow>10n**10n){v.push([val,F(p,p-1n),e,true]);break;}
  v.push([val,val,e,false]);
 }
 return v;
}
if(8n**4n>=7n*5n**4n||(X+2n)*8n**5n>=140n*X*3n**5n)throw Error("fine tail bound");


const out=[];
for(const extra of [67,79,307]){
 const T=[7,13,19,31,extra],R=smooth(T).map(Number).sort((a,b)=>a-b);
 for(const pi of [...T.filter(p=>p%12===1),0]){
 const rows=[...T.map(p=>[p,ints(p,p===pi?"Euler":"ex",false).filter(x=>p!==7||x[2]!==4)]),...R.map(p=>[p,ints(p,"reg",p!==5&&p!==11)])].sort((a,b)=>a[0]-b[0]);
 const n=rows.length,su=Array(n+1);su[n]=F(1);
 for(let i=n-1;i>=0;i--)su[i]=mul(su[i+1],rows[i][1].at(-1)[1]);
 let best=F(0),cross=false;
 function go(i,l,h,es){
  if(cross||cmp(l,F(2))>0||cmp(mul(h,su[i]),F(199,100))<0)return;
  if(i===n){let mp=new Map(rows.map(([p],j)=>[p,es[j]]));if(mp.get(13)===1&&(mp.get(191)||mp.get(653)))return;
  if(cmp(h,best)>0)best=h;if(cmp(h,F(2))>=0)cross=true;return;}
  for(const [l2,h2,e] of rows[i][1])go(i+1,mul(l,l2),mul(h,h2),[...es,e]);
 }
 go(0,F(1),F(1),[]);out.push([extra,R,pi,cross,Number(best[0])/Number(best[1]),best.map(String)]);
 }
}
return out;
}
function emptyTail(extra, cap){

const X=10n**6n,Y2=10n**18n, M=10n**15n, S=[7n,13n,19n,31n,BigInt(extra)];
const F=(a,b=1n)=>[BigInt(a),BigInt(b)],mul=(a,b)=>[a[0]*b[0],a[1]*b[1]],add=(a,b)=>[a[0]*b[1]+a[1]*b[0],a[1]*b[1]],cmp=(a,b)=>a[0]*b[1]-a[1]*b[0],prod=as=>as.reduce(mul,F(1));
function isqrt(n){if(n<2n)return n;let a=1n<<BigInt(Math.ceil(n.toString(2).length/2)),b=(a+n/a)>>1n;while(b<a){a=b;b=(a+n/a)>>1n;}return a;}
function sixth(n){let a=1000n,b=1000n*n;while(b-a>1n){let c=(a+b)>>1n;if(c**6n<=n*10n**18n)a=c;else b=c;}return a;}
const baseProd=prod(S.map(p=>{let z=sixth(p);return F(z,z-1000n);}));
const H=F(...cap);
function tail(q){
 let sum=0n,n=0;
 function walk(i,m){if(i===S.length){let d=isqrt(m);if(d<X)d=X;sum+=M/d+1n;n++;return;}
 for(let w=m;w<=Y2;w*=S[i])walk(i+1,w);}
 for(let w=q;w<=Y2;w*=q)walk(0,w);
 const d=mul(F(X+2n,X),add(F(sum,M),mul(F(1,1000000),mul(baseProd,F(1000n,sixth(q)-1000n)))));
 return [add(d,F(140,10n**8n)),n];
}
function pass(q){let [d]=tail(q);return cmp(add(mul(H,F(q,q-1n)),mul(F(2),d)),F(2))<0;}

let l=522n,h=10000000000n;if(!pass(h))throw Error("bad search");
while(h-l>1n){let m=(h+l)>>1n;if(pass(m))h=m;else l=m;}
let [d,c]=tail(h);return [String(h),Number(d[0])/Number(d[1]),c];

}
const profiles=coreMaxima(), caps=[[67,"9715040803","4867153920",1418],[79,"7465495851947","3733004880000",397858275],[307,"3021371507","1516320000",523]], report=[];
for(const [extra,num,den,expected] of caps){
const applicable=profiles.filter(x=>x[0]===extra);
check(applicable.length===2&&applicable.every(x=>!x[3])&&applicable.some(x=>x[5].join()===num+","+den),"core bound witness");
check(applicable.every(x=>BigInt(x[5][0])*BigInt(den)<=BigInt(num)*BigInt(x[5][1])),"core bound verification");
const ans=emptyTail(extra,[num,den]);
check(+ans[0]===expected,"tail threshold");
report.push([extra,ans]);
}
console.log(JSON.stringify(report));
}

function sevenScreen(suppliedStageInput,console){
// Run with Node.js. Output surviving T for a7=4 and b=2801.
const X=1e6,P=500;
const is=new Uint8Array(X+1).fill(1,2),rem=new Float64Array(X+1),fs=new Map(),ps=[];
for(let p=2;p*p<=X;p++)if(is[p])for(let j=p*p;j<=X;j+=p)is[j]=0;
for(let p=5;p<=X;p++)if(is[p]){ps.push(p);rem[p]=p*p+p+1;fs.set(p,[]);}
function pow(a,e,p){let b=1;while(e){if(e%2)b=b*a%p;a=a*a%p;e=Math.floor(e/2);}return b;}
for(const p of [3,...ps]){let a=[];if(p===3)a=[1];else if(p%3===1){for(let k=2;k<p;k++){let z=pow(k,(p-1)/3,p);if(z!==1){a=[z,p-1-z];break;}}}for(const z of a)for(let j=z;j<=X;j+=p)if(rem[j]){fs.get(j).push(p);while(rem[j]%p===0)rem[j]/=p;}}
const caps=new Map();
for(const p of ps){let f=fs.get(p);if(rem[p]>1)f.push(rem[p]);let l=p%3===1?p:f.at(-1);if(p>P)(caps.get(l)??caps.set(l,[]).get(l)).push(p);}

const F=(a,b=1)=>[BigInt(a),BigInt(b)], mul=(a,b)=>[a[0]*b[0],a[1]*b[1]], add=(a,b)=>[a[0]*b[1]+a[1]*b[0],a[1]*b[1]], cmp=(a,b)=>a[0]*b[1]-a[1]*b[0], prod=as=>as.reduce(mul,F(1)), hp=p=>F(p,p-1);
function req(v,msg){if(!v)throw Error(msg);}
function froot(p){let l=1000,h=p*1000;while(h-l>1){let k=Math.floor((l+h)/2);if(BigInt(k)**6n<=BigInt(p)*10n**18n)l=k;else h=k;}return F(l,l-1000);}

const weightBound=F(1006,1000), W=new Map([...caps].map(([l,a])=>[l,prod(a.map(hp))]));
for(const weight of W.values())req(cmp(weight,weightBound)<0,"tail group below 1.006");
const small=ps.filter(p=>p<=P), other=small.filter(p=>p%3===2&&p!==5&&p!==11);
const cand=[...new Set([...small.filter(p=>p%3===1),...small.filter(p=>p%3===2).flatMap(p=>fs.get(p))])].sort((a,b)=>a-b);

const specs=[
 ["5 exceptional, 11 regular, 31 absent",[7,13,19],1,[31]],
 ["5 exceptional, 11 regular, 31 present",[7,13,19,31],1,[]],
 ["5 regular, 11 exceptional, 19 absent",[7,13,31],1,[19]],
 ["5 regular, 11 exceptional, 19 present",[7,13,31,19],1,[]],
 ["both regular, 13 absent",[7,19,31],0,[13]],
 ["both regular, 13 present",[7,19,31,13],0,[]]];
const candSet=new Set(cand);
const stats=[];
const topWeights=[...W.entries()].sort((a,b)=>{let d=cmp(b[1],a[1]);return d>0n?1:d<0n?-1:0;}).slice(0,12);
for(const [name,fixedT,e,banned] of specs.slice(5)){
 const setup=cand.filter(p=>!fixedT.includes(p)&&!banned.includes(p));const free=6-e-fixedT.length;
 const anonymous=ps.filter(p=>p%3===1&&!candSet.has(p)&&!banned.includes(p)).slice(0,6);
 let n=0,surv=0,best=F(0),worst,remain=[];
 function rec(i,t){
  const T=[...fixedT,...t], S=new Set(T), L=new Set([7,19,31,...T]);
  const reg=other.filter(p=>fs.get(p).every(r=>S.has(r))), opt=other.filter(p=>!reg.includes(p));
  const active=[5,7,11,...reg,...T.filter(p=>p!==7&&p<=P)];
  for(let u=0;u<=free-t.length;u++){
   if(T.length+u<5||(!T.some(p=>p%12===1)&&!u))continue;
   const used=opt.slice(0,free-t.length-u);
   const ub=prod([...active.map(hp),...used.map(hp),...[...L].map(r=>W.get(r)??F(1)),...topWeights.filter(([r])=>!L.has(r)).slice(0,6-L.size).map(([,w])=>w)]);
   const d=add(F(6,X),mul(F(X+2,X*10000),prod([...T,...anonymous.slice(0,u)].map(froot))));
   const threshold=add(ub,mul(F(2),d));
   n++;if(cmp(threshold,F(2))>=0){surv++;remain.push([T,u]);}
   if(cmp(threshold,best)>0){best=threshold;worst=[T,u,active,used,Number(ub[0])/Number(ub[1]),Number(d[0])/Number(d[1])];}
  }
  if(t.length<free)for(let j=i;j<setup.length;j++)rec(j+1,[...t,setup[j]]);
 }
 rec(0,[]);
 stats.push({name,n,surv,threshold:Number(best[0])/Number(best[1]),worst,remain});
} 




const choices=stats[0].remain.filter(([T,u])=>T.length===5&&!u);
const keep=[];
const hp7=F(2801,2401), b=2801;
for(const [T] of choices){
 const S=new Set(T), R=ps.filter(p=>p%3===2&&p!==b&&fs.get(p).every(r=>S.has(r)));
 const U=prod([...T.filter(p=>p!==7),...R,b].map(hp).concat([hp7]));
 const d=mul(F(X+2,X*10000),prod(T.map(froot)));
 if(cmp(add(U,mul(F(2),d)),F(2))>=0n)keep.push([T,R]);
}
if(choices.length!==73||keep.map(([T])=>T.at(-1)).join(',')!=='67,79,307')throw Error('a7=4 support screen');
console.log(JSON.stringify(keep.map(([T])=>T)));
}

function sevenIntervals(suppliedStageInput,console){
const X=10n**16n, max=X*X+X+1n;
const G=([a,b])=>{let x=a,y=b;while(y)[x,y]=[y,x%y];return [a/x,b/x];};
const F=(a,b=1n)=>[BigInt(a),BigInt(b)], mul=(a,b)=>G([a[0]*b[0],a[1]*b[1]]), cmp=(a,b)=>a[0]*b[1]-a[1]*b[0];
function isqrt(n){if(n<2n)return n;let g=1n<<BigInt(Math.ceil(n.toString(2).length/2));let h=(g+n/g)>>1n;while(h<g){g=h;h=(g+n/g)>>1n;}return g;}
const ok=Array.from({length:576},(_,i)=>i).reduce((s,i)=>(s.add(i*i%576),s),new Set());
function smooth(T){
 let ans=[];
 function go(i,v){
  if(i===T.length){let d=4n*v-3n;if(!ok.has(Number(d%576n)))return;let r=isqrt(d),q=(r-1n)/2n;if(r*r===d&&q%6n===5n&&q<=X)ans.push(q);return;}
  for(let w=v;w<=max;w*=BigInt(T[i]))go(i+1,w);
 }
 go(0,1n);return ans;
}
let exponent7=0;
function ints(p,role,optional){
 p=BigInt(p);if(p===7n){let den=7n**BigInt(exponent7),v=F((7n*den-1n)/6n,den);return [[v,v]];}const v=optional?[[F(1),F(1)]]:[];
 for(let e=role==="Euler"?1:2;;e+=role==="Euler"?4:2){
  if(((e+1)%3===0)!==(role==="reg"))continue;
  let pow=p**BigInt(e),val=F(p*pow-1n,(p-1n)*pow);
  if((p-1n)*pow>10n**10n){v.push([val,F(p,p-1n)]);break;}
  v.push([val,val]);
 }
 return v;
}
if(8n**4n>=7n*5n**4n||(X+2n)*8n**5n>=140n*X*3n**5n)throw Error("fine tail bound");
const hiThreshold=F(2); // 2(1 - 140/10^8)

function req(condition,message){if(!condition)throw Error(message);}


function order(p,R){
 const r=BigInt(R), x=p%r; if(!x)throw Error("nonunit");
 const pow=(e)=>{let a=x,b=1n;for(let n=BigInt(e);n;n>>=1n,a=a*a%r)if(n&1n)b=b*a%r;return b};
 let ord=R-1,k=ord;
 for(let d=2;d*d<=k;d++)if(k%d===0){while(k%d===0)k/=d;while(ord%d===0&&pow(ord/d)===1n)ord/=d;}
 if(k>1)while(ord%k===0&&pow(ord/k)===1n)ord/=k;
 return ord;
}

let forcedRegular=new Set();
function regular(T,p,optional){
 const es=[...new Set(T.map(r=>order(p,r)).filter(n=>n%6===3))].sort((a,b)=>a-b).map(n=>n-1);
 const choices=optional&&!forcedRegular.has(p)?[[F(1),F(1)]]:[];
 let cap=0,maxpow=1n;while((p-1n)*maxpow<=10n**10n){cap++;maxpow*=p;}
 for(const e of es){
  let power=p**BigInt(Math.min(e,cap)), h=F(p*power-1n,(p-1n)*power);
  if(e>=cap){choices.push([h,F(p,p-1n)]);break;}
  choices.push([h,h]);
 }
 return choices;
}
function prime(p){for(let j=2;j*j<=p;j++)if(p%j===0)return false;return true;}
function boxes(T,pi,R,q,B){
 const count=T.length;
 const fine=count===5?140n:360n;
 if((X+2n)*8n**BigInt(count)>=fine*X*3n**BigInt(count))throw Error("fine tail bound");
 const loThreshold=F(2n*(10n**8n-fine),10n**8n);
 const domains=[
  ...T.map(p=>ints(p,p===pi?"Euler":"ex",false)),
  ...(q?[ints(q,"ex",false)]:[]),
  ...R.filter(p=>BigInt(q)!==p).map(p=>regular(T,p,p!==5n&&p!==11n))
 ];
 if(domains.some(x=>!x.length))return null;
 const suffix=Array(domains.length+1);suffix[domains.length]=q||T.length===6?F(1):F(B+1,B);
 for(let i=domains.length-1;i>=0;i--)suffix[i]=mul(suffix[i+1],domains[i].at(-1)[1]);
 let best=F(0),num=0;
 function go(i,lo,hi){if(num)return;
  if(cmp(lo,hiThreshold)>0||cmp(mul(hi,suffix[i]),loThreshold)<0)return;
  if(i===domains.length){num++;if(cmp(hi,best)>0)best=hi;return;}
  for(const [a,b] of domains[i])go(i+1,mul(lo,a),mul(hi,b));
 }
 go(0,F(1),F(1));if(!num)return null;
 const gap=loThreshold[0]*best[1]-loThreshold[1]*best[0];
 return [q,num,gap>0n?String((loThreshold[0]*best[1])/gap):"unbounded"];
}




function check(ok,why){if(!ok)throw Error(why);}
const sig=(p,a)=>(p**BigInt(a+1)-1n)/(p-1n);
const phi=p=>p*p+p+1n;
const primeCheck=nums=>nums.every(p=>prime(Number(p)));
check(sig(7n,4)===2801n&&phi(2801n)===37n*43n*4933n&&primeCheck([2801,37,43,4933]),"a7=4 factors");
check(sig(7n,6)===29n*4733n&&phi(29n)===13n*67n&&phi(4733n)===22406023n&&primeCheck([29,4733,67,22406023]),"a7=6 factors");
check(sig(7n,10)===1123n*293459n&&phi(293459n)===277n*310897033n&&primeCheck([1123,293459,277,310897033]),"a7=10 factors");
check(sig(7n,12)===16148168401n &&sig(7n,16)===14009n*2767631689n&&primeCheck([16148168401,14009,2767631689]),"a7=12,16 factors");
check(16148168401%3===1&&2767631689%3===1,"large forced targets");
for(const r of [16148168401,2767631689])for(let p=5;p<=1000000;p++)check((p*p+p+1)%r!==0,"large forced target fell in a nonempty group");
const big=4534166740403n;
check(sig(7n,18)===419n*big&&phi(419n)===13n*13537n&&primeCheck([419,big,13537,199,229,373]),"a7=18 factors");
check(phi(big)%(199n*229n*373n)===0n&&[199,229,373].every(p=>p%3===1),"too many targets from big regular base");
const alternatives=readInput().map(T=>[4,T,2801,[]]);
alternatives.push([6,[7,13,19,31,67],4733,[29]],[6,[7,13,19,31,67,22406023],0,[29,4733]],[10,[7,13,19,31,1123],293459,[]],[18,[7,13,19,31,13537],Number(big),[419]]);
let count=0;
for(const [a,T,b,mandatory] of alternatives){
 exponent7=a;
 const R=smooth(T).filter(p=>T.every(t=>p%BigInt(t)!==0n));
 forcedRegular=new Set(mandatory.map(BigInt));
 if([...forcedRegular,5n,11n].some(p=>!R.includes(p)))throw Error("missing regular base");
 for(const pi of T.filter(p=>p%12===1)){
   count++;check(!boxes(T,pi,R,b,500),"surviving Euler choice "+a+" "+pi);
 }
}
check(count===8,"Euler choices");
console.log(JSON.stringify({caseCount:count,eliminated:[4,6,10,12,16,18]}));

function readInput(){return suppliedStageInput;}
}

const noThirteen=stage(thirteenIntervals,stage(thirteenScreen));
const anonymous=stage(anonymousScreen);
const noTwentyThree=stage(excludeTwentyThree);
const noSeventeen=stage(seventeenIntervals,stage(seventeenScreen));
const six=stage(sixThresholds);
const seven=stage(sevenIntervals,stage(sevenScreen));
console.log(JSON.stringify({noThirteen,anonymous:anonymous.map(x=>[x.T,x.remaining]),noTwentyThree,
 noSeventeen:{eulerChoices:noSeventeen.eulerChoices,excluded:noSeventeen.excluded},six,seven},null,2));
