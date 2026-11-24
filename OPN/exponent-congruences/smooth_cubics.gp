\\ Complete positive q = 5 mod 6 with q^2+q+1 supported on the given primes.
\\ Run with PARI/GP: gp -q -f smooth_cubics.gp

eis_mul(z,w) = [z[1]*w[1]-z[2]*w[2], z[1]*w[2]+z[2]*w[1]-z[2]*w[2]];

split_choices(p) = {
  my(bound=sqrtint(4*p)+1, pi, pc);
  if(!isprime(p) || p%3 != 1, error("expected a split rational prime"));
  forstep(b=3,bound,3,
    for(a=-bound,bound,
      if((a-1)%3 == 0 && a^2-a*b+b^2 == p,
        pi=[a,b]; pc=[a-b,-b];
        return([[1,0], pi, pc, eis_mul(pi,pi), eis_mul(pc,pc)])
      )
    )
  );
  error("primary factor missing");
};

smooth_bases(targets) = {
  my(n=#targets, options=vector(n,i,split_choices(targets[i])));
  my(after2=0, equations=0, hits=List());
  forvec(v=vector(n,i,[1,5]),
    my(gamma=[-1,-1], g, h, valid=1, poly, sols);
    for(i=1,n, gamma=eis_mul(gamma,options[i][v[i]]));
    g=gamma[1]; h=gamma[2];
    if(g%2 != 1 || h%2 != 1, next);
    after2++;
    for(i=1,n,
      if(v[i] != 1 && Mod(h,targets[i])^((targets[i]-1)/3) != 1, valid=0)
    );
    if(!valid, next);
    equations++;
    poly=h*'x^3+3*(g-h)*'x^2-3*g*'x+h;
    sols=thue(thueinit(poly,1),-1);
    for(i=1,#sols,
      my(z=sols[i], value=eis_mul(gamma,eis_mul(z,eis_mul(z,z))));
      my(q=value[1], rem=q^2+q+1);
      if(value[2] != -1, error("incorrect Thue solution"));
      if(q<=0 || q%6 != 5, next);
      for(j=1,n, rem=rem/targets[j]^valuation(rem,targets[j]));
      if(rem==1, listput(hits,q))
    )
  );
  return([after2,equations,Set(Vec(hits))]);
};

verify_tables() = {
  my(base=[7,13,19,31], extras=[0,67,79,307]);
  my(expected=[[205,9,[5,11,191,653]], [1045,24,[5,11,29,191,653]],
               [1045,19,[5,11,23,191,653]], [1025,25,[5,11,17,191,653]]]);
  for(i=1,#extras,
    my(targets=if(extras[i],concat(base,extras[i]),base));
    my(result=smooth_bases(targets));
    if(result != expected[i], error("table mismatch for ",targets));
    print(targets,"  ",result)
  );
};

verify_tables();
