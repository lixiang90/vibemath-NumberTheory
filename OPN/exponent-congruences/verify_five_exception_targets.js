// Node.js standard library only. Exact five non-Euler exceptions with 3 absent.
let data;
{
const X=1e6,P=200;
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
for(const [name,fixedT,e,banned] of specs){
 const setup=cand.filter(p=>!fixedT.includes(p)&&!banned.includes(p));const free=6-e-fixedT.length;
 const anonymous=ps.filter(p=>p%3===1&&!candSet.has(p)&&!banned.includes(p)).slice(0,6);
 let n=0,surv=0,best=F(0),worst,remain=[];
 function rec(i,t){
  const T=[...fixedT,...t], S=new Set(T), L=new Set([7,19,31,...T]);
  const reg=other.filter(p=>fs.get(p).every(r=>S.has(r))), opt=other.filter(p=>!reg.includes(p));
  const active=[5,7,11,...reg,...T.filter(p=>p!==7&&p<=P)];
  for(let u=0;u<=free-t.length;u++){
   if(!T.some(p=>p%12===1)&&!u)continue;
   const used=opt.slice(0,free-t.length-u);
   const ub=prod([...active.map(hp),...used.map(hp),...[...L].map(r=>W.get(r)??F(1)),...Array(6-L.size).fill(weightBound)]);
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


req(cand.length===35,"explicit target set");
const four=stats.map(x=>x.remain.filter(([T,u])=>T.length+u===4));
req(four.every(rows=>rows.every(([,u])=>!u)),"anonymous four-target support");
req(four.map(rows=>rows.length).join(",")==="12,1,3,1,7,1","four-target screen counts");
const five=stats.slice(0,4).map((s,k)=>s.remain.filter(([T,u])=>{
 if(T.length+u!==5)return false;
 req(!u,"anonymous five-target support");
 const forced=k<2?5:11,S=new Set(T);
 const regular=ps.filter(p=>p%3===2&&p!==forced&&fs.get(p).every(q=>S.has(q)));
 const U=prod([...T,forced,...regular].map(hp));
 const D=add(F(6,X),mul(F(X+2,X*10000),prod(T.map(froot))));
 return cmp(add(U,mul(F(2),D)),F(2))>=0;
}).map(([T])=>T));
req(five.map(rows=>rows.length).join(",")==="19,3,3,3","five-target screen counts");

// Exactly three targets and a single exception among 5 and 11.
const pairs=[];
for(const [T,forced,required] of [[[7,13,19],5,11],[[7,13,31],11,5]]){
 const S=new Set(T);
 const regular=ps.filter(p=>p%3===2&&p!==forced&&fs.get(p).every(t=>S.has(t)));
 const optional=ps.filter(p=>p%3===2&&p!==forced&&p!==required);
 // Double counting optional exceptions here only enlarges the cap.
 const U=prod([...T,forced,...regular].map(hp));
 const D=add(F(6,X),mul(F(X+2,X*10000),prod(T.map(froot))));
 const keep=x=>cmp(add(mul(U,x),mul(F(2),D)),F(2))>=0;
 req(!keep(hp(optional[0])),"cannot yet ignore unknown exceptions above X");
 for(let i=0;i+1<optional.length;i++){
  const a=optional[i];
  if(!keep(mul(hp(a),hp(optional[i+1]))))break;
  for(let j=i+1;j<optional.length;j++){
   const b=optional[j];
   if(!keep(mul(hp(a),hp(b))))break;
   pairs.push([T,[forced,a,b],[required]]);
  }
 }
}
req(pairs.length===4,"three-target exception pairs");
data={four:four.map(rows=>rows.map(([T])=>T)),five,pairs,smallExceptions:ps.filter(p=>p<=P&&p%3===2)};
}
{
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
const hiThreshold=F(2),loThreshold=F(199999720n,100000000n); // 2(1 - 140/10^8)

function req(condition,message){if(!condition)throw Error(message);}

const smoothCache=new Map();
function boxes(T,pi,exceptions,required,unknown){
 const key=[...T].sort((a,b)=>a-b).join(',');
 if(!smoothCache.has(key))smoothCache.set(key,smooth(T));
 const regular=smoothCache.get(key).filter(p=>!exceptions.some(q=>p===BigInt(q)));
 req(required.every(p=>regular.includes(BigInt(p))),"required regular preimage");
 const domains=[
  ...T.map(p=>ints(p,p===pi?"Euler":"ex",false)),
  ...exceptions.map(p=>ints(p,"ex",false)),
  ...regular.map(p=>ints(p,"reg",!required.some(q=>p===BigInt(q))))
 ];
 const suffix=Array(domains.length+1);
 // Every unlocated exceptional base exceeds 200 and contributes < 201/200.
 suffix[domains.length]=F(201n**BigInt(unknown),200n**BigInt(unknown));
 for(let i=domains.length-1;i>=0;i--)suffix[i]=mul(suffix[i+1],domains[i].at(-1)[1]);
 function go(i,lo,hi){
  if(cmp(lo,hiThreshold)>0||cmp(mul(hi,suffix[i]),loThreshold)<0)return false;
  return i===domains.length||domains[i].some(([a,b])=>go(i+1,mul(lo,a),mul(hi,b)));
 }
 return go(0,F(1),F(1));
}

for(const [T,exceptions,required] of data.pairs){
 for(const pi of T.filter(q=>q%12===1))req(!boxes(T,pi,exceptions,required,0),"three-target case survived");
}
let Euler5=0;
for(let k=0;k<4;k++){
 const exceptions=[k<2?5:11],required=[k<2?11:5];
 for(const T of data.five[k])for(const pi of T.filter(q=>q%12===1)){
  Euler5++;
  req(!boxes(T,pi,exceptions,required,0),"five-target case survived");
 }
}
req(Euler5===37,"five-target Euler choices");

const remain=[];
let Euler4=0;
for(let k=0;k<6;k++){
 const fixed=k<2?[5]:k<4?[11]:[];
 const required=k<2?[11]:k<4?[5]:[5,11];
 const optional=data.smallExceptions.filter(p=>!fixed.includes(p)&&!required.includes(p));
 const slots=2-fixed.length;
 for(const T of data.four[k])for(const pi of T.filter(q=>q%12===1)){
  Euler4++;
  function enumerate(i,known){
   const count=slots-known.length;
   if(boxes(T,pi,[...fixed,...known],required,count))remain.push({k,T:[...T].sort((a,b)=>a-b),pi,known,count});
   if(count)for(let j=i;j<optional.length;j++)enumerate(j+1,[...known,optional[j]]);
  }
  enumerate(0,[]);
 }
}
req(Euler4===30,"four-target Euler choices");
req(remain.length===2&&remain.every((r,i)=>r.k===5&&r.T.join(",")==="7,13,19,31"&&
 r.pi===13&&r.count===1&&r.known.join(",")===(i===0?"17":"23")),
 "unexpected four-target survivors");
console.log(JSON.stringify({threeTargetPairs:data.pairs,fiveTargetEulerChoices:Euler5,
 fourTargetEulerChoices:Euler4,remainingFourTargetCases:remain},null,2));
}
