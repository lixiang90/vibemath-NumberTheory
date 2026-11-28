// Run with Node.js. The merged PR #2 supplies the initial four-target reduction.
function screen(cutoff){
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
const hiThreshold=F(2),loThreshold=F(199999720n,100000000n); // 2(1 - 140/10^8)

const T=[7,13,19,31], reg=smooth(T).map(Number).sort((a,b)=>a-b);
const B=Number(cutoff), stats=[];
function need(ok,message){if(!ok)throw Error(message);}
need(reg.join(",")==="5,11,191,653","regular bases through Y");
function regularChoices(p){let base=BigInt(p),c=F(base*base+base+1n,base*base);return p===5||p===11?[[c,c,2,false]]:[[F(1),F(1),0,false],[c,c,2,false]];}
function sieve(n){const a=new Uint8Array(n+1).fill(1,2);for(let j=2;j*j<=n;j++)if(a[j])for(let k=j*j;k<=n;k+=j)a[k]=0;return [...a.keys()].filter(k=>a[k]&&k%3===2&&k>200);}
for(const fixed of [17,23])for(const q of [...sieve(B),0]){
 const E=q?[fixed,q]:[fixed];
 const choices=[...T.map(p=>[p,ints(p,p===13?"Euler":"ex",false)]),...E.map(p=>[p,ints(p,"ex",false)]),...reg.filter(p=>!E.includes(p)).map(p=>[p,regularChoices(p)])].sort((a,b)=>a[0]-b[0]);
 const suffix=Array(choices.length+1);suffix[choices.length]=q?F(1):F(B+1,B);
 for(let i=choices.length-1;i>=0;i--)suffix[i]=mul(suffix[i+1],choices[i][1].at(-1)[1]);
 const remain=[];let best=F(0);
 function go(i,l,h,es){if(cmp(l,hiThreshold)>0||cmp(mul(h,suffix[i]),loThreshold)<0)return;
 if(i===choices.length){const map=new Map(es.map(([p,e,s])=>[p,[e,s]]));if(map.get(13)[0]===1&&[191,653].some(p=>map.get(p)[0]>0))return;if(cmp(h,best)>0)best=h;remain.push(es);return;}
 for(const [a,b,e,sat] of choices[i][1])go(i+1,mul(l,a),mul(h,b),[...es,[choices[i][0],e,sat]]);
 }
 go(0,F(1),F(1),[]);
 if(remain.length)stats.push([fixed,q,remain,best]);
}


return stats;
}

function need(ok,message){if(!ok)throw Error(message);}
function ord(base,p){
 let value=1;
 for(let k=1;k<p;k++){value=value*base%p;if(value===1)return k;}
 throw Error("undefined order");
}
const targets=[7,13,19,31], bases=[5,11,191,653];
const allowed=bases.map(p=>targets.map(r=>ord(p,r)).filter(order=>order%3===0&&order%2===1).sort((a,b)=>a-b).join(","));
need([...new Set(targets.flatMap(r=>Array.from({length:r-1},(_,i)=>i+1).filter(n=>n%6===3&&(r-1)%n===0)))].sort((a,b)=>a-b).join(',')==='3,9,15',"possible length count");
need(allowed.join(";")==="3,9;3,3;3,3,3;3,3,3","check regular exponent orders");
need(5**6+5**3+1===19*829,"Phi9(5)");
for(let j=2;j*j<=829;j++)need(829%j!==0,"829 prime");
for(const p of [191,653])need((p*p+p+1)%169===0,"13-adic obstruction");
const first=screen(3000);
need(first.length===1&&first[0][0]===17&&first[0][1]===0&&first[0][2].length===72,"first screen");
const [Hn,Hd]=first[0][3], Ln=199999720n, Ld=100000000n;
const gap=Ln*Hd-Ld*Hn;
need(gap>0n&&(Ln*Hd)/gap===3729n,"upper bound on second exceptional prime");
const second=screen(3729);
need(second.map(([p,q,boxes])=>p+","+q+","+boxes.length).join(";")==="17,3041,18;17,3701,28;17,3719,28","second screen");
need(ord(13,17)===4,"Euler contributor test");
for(const q of [7,19,31,3041,3701,3719])need(ord(q,17)%2===0,"even-powered exceptional contributor test");
console.log(JSON.stringify({maximumKnownAbundance:[String(Hn),String(Hd)],boundOnSecondExceptionalPrime:3729,discarded:second.map(([p,q,boxes])=>[p,q,boxes.length]),result:"four-target branch excluded"}));

