def approx(expr,N):
    """
    Заменяет oo на N в любом символьном выражении
    """
    import re
    return SR(re.sub("Infinity", str(N), str(expr))).unhold()

class TrigonometricSeries:
    
    def __init__(self,expr):
        self.expr=expr
        [term,n,n_min,n_max]=self.expr.operands()
        self.term = term
        self.n=n
        self.min=n_min
        self.max=n_max
        L=self.list()
        try:
            [a,b,x]=self.list()
            self.a=a
            self.b=b
            self.x=x
        except:
            return None
        
    def latex(self):
        return latex(self.expr)
    
    def list(self):
        if str(self.expr.operator())!='sum':
            return 'Error: the expression is not a sum'
        term = self.term
        n = self.n
        a = SR.wild(0)
        b = SR.wild(1) 
        x = SR.wild(2)
        D=term.match(a*cos(n*x)+b*sin(n*x))
        if D is not None:
            for xx in D[x].variables():
                if diff(D[a],xx)!=0 or diff(D[b],xx)!=0:
                    return 'Error: the expression is not a trigonometric series'
            return [D[a],D[b],D[x]]
        D=term.factor().match(a*cos(n*x))
        if D is not None:
            for xx in D[x].variables():
                if diff(D[a],xx)!=0:
                    return 'Error: the expression is not a trigonometric series'
            return [D[a],0,D[x]]
        D=term.factor().match(b*sin(n*x))
        if D is not None:
            for xx in D[x].variables():
                if diff(D[b],xx)!=0:
                    return 'Error: the expression is not a trigonometric series'
            return [0,D[b],D[x]]
        else:
            return 'Error: the expression is not a trigonometric series'
        
    def kryloff_slow_part(self,r,list=False):
        a_slow=taylor(self.a,n,oo,r)
        b_slow=taylor(self.b,n,oo,r)
        if list==False:
            return sum(a_slow*cos(self.n*x)+b_slow*sin(self.n*x),self.n,self.min,self.max, hold=True)
        else:
            return [a_slow, b_slow]
      
    def eval(self, method='bernoulli'):
        if method=='bernoulli':
            if self.min!=1 or self.max!=oo:
                return 'Error: the limits in the sum must be 1 and oo'
            f=0
            if self.a!=0:
                for [c,p] in (self.a).coefficients():
                    if diff(c,self.n)!=0:
                        return None
                    elif p==-1:
                        f=f-c*ln(2*abs(sin(x/2)))
                    elif p in ZZ and -p>0 and ZZ(p) % 2 ==0:
                        f=f+c*self.bernoulli_p(x,-p)
                    else:
                        f=f+c*function('kryloff_cos')(x,-p)
            if self.b!=0:
                for [c,p] in (self.b).coefficients():
                    if diff(c,self.n)!=0:
                        return None
                    elif p in ZZ and -p>0 and ZZ(p) % 2 ==1:
                        f=f+c*self.bernoulli_p(x,-p)
                    else:
                        f=f+c*function('kryloff_sin')(x,-p)   
            return f
        else:
            return sum(self.term,self.n,self.min,self.max,hold=False)
    
    def bernoulli_p(self,x,n):
        b=bernoulli_polynomial(x/2/pi,n)
        d=(-1)^(-floor(n/2)-1)*(2*pi)^n*b/2/factorial(n)
        return SR(d)
    
    def diff(self,u):
        [a_slow,b_slow]=self.kryloff_slow_part(1, list=True)
        slow_part=sum(a_slow*cos(self.n*x)+b_slow*sin(self.n*x),self.n,self.min,self.max, hold=True)
        g=((self.a-a_slow)*cos(self.n*x).diff(u)).factor() + ((self.b-b_slow)*sin(self.n*x).diff(u)).factor()
        if g==0:
            fast_part_diff=0
        else:
            fast_part_diff=sum(g,self.n,self.min,self.max, hold=True)
        return TrigonometricSeries(slow_part).eval().diff(u) + fast_part_diff
    
    def ode(self, field=QQ):
        K=FractionField(field[self.n])
        D1=SR(K(self.b).denominator())
        D2=D1.subs(self.n==-self.n)
        N1=SR(K(self.b).numerator())
        L=(D1*D2).expand()
        bb=(N1*D2).expand()
        var('D')
        return [L.subs(n==i*D), sum(bb*sin(self.n*self.x),self.n,self.min,self.max, hold=True)]
        
