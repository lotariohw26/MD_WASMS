import sympy, pandas, numpy, math, polysolver
import pandas as pd
import polysolver as plsv
import math                                                                   
from sympy import solve, Eq, symbols, latex, simplify, diff, poly, sympify, Matrix, pprint, collect, expand, Poly, Symbol, Pow
def rall(sel=[0, 0, 0]):
    n, m = symbols('n m')
    
    rxy = Matrix([[m, -n, 0], 
                  [n, m, 0],
                  [0, 0, 1]])
    
    rxz = Matrix([[m, 0, -n], 
                  [0, 1, 0],
                  [n, 0, m]])
    
    ryx = Matrix([[m, n, 0], 
                  [-n, m, 0],
                  [0, 0, 1]])
    
    ryz = Matrix([[1, 0, 0], 
                  [0, m, -n],
                  [0, n, m]])
    
    rzx = Matrix([[m, 0, n], 
                  [0, 1, 0],
                  [-n, 0, m]])
    
    rzy = Matrix([[1, 0, 0], 
                  [0, m, n],
                  [0, -n, m]])
    
    ps = [0, rxy, rxz, ryx, ryz, rzx, rzy]
    
    allrot = [ps[i] for i in sel]
    return allrot

def genpolycoeffn(form,expr,solv):
    forp = {'N': 1, 'H': 2, 'O': 3}
    elem = [["alpha","x","y","zeta","lamda","Omega"], ["alpha","g","h","Gamma","Omega","lamda"], ["alpha","m","n","xi","lamda","Omega"]]
    x, y, z = sympy.symbols('x y z')
    alpha, g, h, n, m, zeta, Gamma, lamda, ui = symbols('alpha g h n m zeta Gamma lamda ui')
    k0, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 = symbols('k0:11')
    parm = elem[2-1][:5]
    ls, rs = expr.split('=')
    expr = Eq(sympify(ls), sympify(rs))
    plr = int(sympy.total_degree(expr))-1
    polys = poly(expr.rhs - expr.lhs, sympify(solv)).all_coeffs()
    abc = [1,0,0,0,1,0,0,0,1]
    uvw = []
    ABCDE = [0, 0, 0, 0, 0]
    ABCDE[:len(polys)] = polys
    matarch = 0
    return ABCDE, abc, matarch, plr

#genpolycoeffn(form=2,expr="alpha=k0+k1*g+k2*h",solv='g')

def genpolycoeffr(elem,expr,solv,eur):
    #elem = [["alpha","x","y","zeta","lamda","Omega"], ["alpha","g","h","Gamma","Omega","lamda"], ["alpha","m","n","xi","lamda","Omega"]]
    #breakpoint()
    x, y, z = sympy.symbols('x y z')
    alpha, g, h, n, m, zeta, Gamma, lamda, ui = symbols('alpha g h n m zeta Gamma lamda ui')
    k0, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 = symbols('k0:11')
    parm = elem #elem[form-1][:5]
    pamm = [element + "_m" for element in parm]
    ls, rs = expr.split('=')
    expr = Eq(sympify(ls), sympify(rs))
    plr = int(sympy.total_degree(expr))-1
    dxyz = {'x': 1, 'y': 2, 'z': 3}
    nrs = parm.index(solv)
    # Defining
    x,  y,  z  = sympy.symbols('x y z')
    ui, vi, wi = sympy.symbols('ui vi wi')
    u0, v0, w0 = sympy.symbols('u0 v0 w0')
    u1, v1, w1 = sympy.symbols('u1 v1 w1')
    u2, v2, w2 = sympy.symbols('u2 v2 w2')
    a1, a2, a3 = sympy.symbols('a1 a2 a3')
    b1, b2, b3 = sympy.symbols('b1 b2 b3')
    c1, c2, c3 = sympy.symbols('c1 c2 c3')
    m1, m2, m3 = sympy.symbols('m1 m2 m3')
    n1, n2, n3 = sympy.symbols('n1 n2 n3')
    mu, mv, mw = sympy.symbols('mu mv mw')
    # Replace
    PS = rall(sel=eur)
    r0 = Matrix([u0, v0, w0])
    r1 = PS[0].subs([(n, n1), (m, m1)]) * r0
    r2 = PS[1].subs([(n, n2), (m, m2)]) * r1
    RR = PS[2].subs([(n, n3), (m, m3)]) * r2
    sympy.expand(RR.row(0)[0])
    sympy.expand(RR.row(1)[0])
    sympy.expand(RR.row(2)[0])
    av = sympy.collect(sympy.expand(RR.row(0)[0]), (u0, v0, w0))
    bv = sympy.collect(sympy.expand(RR.row(1)[0]), (u0, v0, w0))
    cv = sympy.collect(sympy.expand(RR.row(2)[0]), (u0, v0, w0))
    a1s = av.coeff(u0)
    a2s = av.coeff(v0)
    a3s = av.coeff(w0)
    b1s = bv.coeff(u0)
    b2s = bv.coeff(v0)
    b3s = bv.coeff(w0)
    c1s = cv.coeff(u0)
    c2s = cv.coeff(v0)
    c3s = cv.coeff(w0)
    Eu = Eq(x, a1s * u0 + a2s * v0 + a3s * w0)
    Ev = Eq(y, b1s * u0 + b2s * v0 + b3s * w0)
    Ew = Eq(z, c1s * u0 + c2s * v0 + c3s * w0)
    abc = [a1s, a2s, a3s, b1s, b2s, b3s, c1s, c2s, c3s]
    # Manual
    pvar = [w0, v0, u0, w0**2, v0*w0, v0**2, u0*w0, u0*v0, u0**2, w0**3, v0*w0**2, v0**2*w0, v0**2*w0, v0**3, u0*w0**2, u0*v0*w0, u0*v0**2, u0**2*w0, u0**2*v0, u0**3]
    # Connect to rotmat
    Eu = Eq(x, a1 * u0 + a2 * v0 + a3 * w0)
    Ev = Eq(y, b1 * u0 + b2 * v0 + b3 * w0)
    Ew = Eq(z, c1 * u0 + c2 * v0 + c3 * w0)
    exprf = solve(expr,z)[0]-z
    exprc = exprf.subs([(x, solve(Eu, x)[0]), (y, solve(Ev, y)[0]), (z, solve(Ew, z)[0])])
    # Reorganize
    exprr = collect(expand(exprc), pvar)
    xr = [] 
    yr = []
    zr = []
    clma = ['u0','v0','w0','d','var','expr','expr2']
    eqsn = len(exprr.args)
    data = [["0"] * 7 for _ in range(eqsn)]
    matarch = pandas.DataFrame(data, columns=clma)
    matarch.loc[0, ['expr','expr2']]=exprr.args[0]
    matarch.loc[0, ['d']]='d_000'
    for i in range(1, eqsn):
        expt = exprr.args[i]
        expr = expt.args[-1]
        varn = expt / expr
        nrfs = len(varn.free_symbols)
        #varn.free_symbols
        matarch.loc[i, 'var'] = str(varn)
        for j in range(0,nrfs):
            varn = expt.args[j].as_base_exp()[0]
            pown = expt.args[j].as_base_exp()[1]
            matarch.loc[i, str(varn)] = pown
            subd = {a1: a1s, a2: a2s, a3: a3s,b1: b1s, b2: b2s, b3: b3s,c1: c1s, c2: c2s, c3: c3s}
            matarch.loc[i, 'expr'] = expr
            matarch.loc[i, 'expr2'] = expr.subs(subd) 
            cmbx = int(matarch.loc[i, 'u0'])
            cmby = int(matarch.loc[i, 'v0'])
            cmbz = int(matarch.loc[i, 'w0'])
            dr = sum([cmbx, cmby, cmbz])
            matarch.loc[i,'d']='d_'+"".join([str(dr),str(cmbx),str(cmby)])
            dic1 = matarch.set_index('d')['expr'].to_dict()
            dic2 = matarch.set_index('d')['expr2'].to_dict()
            dic = [dic1,dic2][0]
    ABCDE = [0, 0, 0, 0, 0]
    A=Matrix([0,0,0])
    B=Matrix([0,0,0])
    C=Matrix([0,0,0])
    D=Matrix([0,0,0])
    E=Matrix([0,0,0])
    if plr in [3]:
        A[0] += dic['d_330']
        A[1] += dic['d_303']
        A[2] += dic['d_300']
        B[0] += dic['d_320']*z + dic['d_321']*y
        B[1] += dic['d_302']*z + dic['d_312']*x
        B[2] += dic['d_301']*y + dic['d_310']*x
        C[0] += dic['d_310']*z**2 + dic['d_311']*y*z + dic['d_312']*y**2
        C[1] += dic['d_301']*z**2 + dic['d_311']*x*z + dic['d_321']*x**2
        C[2] += dic['d_302']*y**2 + dic['d_311']*x*y + dic['d_320']*x**2
        D[0] += dic['d_300']*y**3 + dic['d_301']*y*z**2 + dic['d_302']*y**2*y + dic['d_303']*z**3
        D[1] += dic['d_300']*x**3 + dic['d_310']*x*y**2 + dic['d_320']*y**2*x + dic['d_330']*y**3
        D[2] += dic['d_303']*y**3 + dic['d_312']*z*y**2 + dic['d_321']*z**2*y + dic['d_330']*y**3
    if plr in [2,3]:
        B[0] += dic['d_220']
        B[1] += dic['d_202']
        B[2] += dic['d_200']
        C[0] += dic['d_210']*z + dic['d_211']*y
        C[1] += dic['d_201']*z + dic['d_211']*x
        C[2] += dic['d_201']*y + dic['d_210']*x
        D[0] += dic['d_200']*z**2 + dic['d_201']*y*z + dic['d_202']*y**2
        D[1] += dic['d_200']*z**2 + dic['d_210']*x*z + dic['d_220']*x**2
        D[2] += dic['d_202']*y**2 + dic['d_211']*x*y + dic['d_220']*x**2
    if plr in [1,2,3]:
        C[0] += dic['d_110']
        C[1] += dic['d_101']
        C[2] += dic['d_100']
        D[0] += dic['d_100']*z + dic['d_101']*y + dic['d_000']
        D[1] += dic['d_100']*z + dic['d_110']*x + dic['d_000']
        D[2] += dic['d_101']*y + dic['d_110']*x + dic['d_000']
    ABCDE[0] = A[nrs].subs([(x,pamm[0]),(y,pamm[1]),(z,pamm[2])])
    ABCDE[1] = B[nrs].subs([(x,pamm[0]),(y,pamm[1]),(z,pamm[2])])
    ABCDE[2] = C[nrs].subs([(x,pamm[0]),(y,pamm[1]),(z,pamm[2])])
    ABCDE[3] = D[nrs].subs([(x,pamm[0]),(y,pamm[1]),(z,pamm[2])])
    ABCDE[4] = E[nrs]
    indpr = [[2,3,4,0,1],[1,2,3,4,0],[0,1,2,3,4]]
    ABCDE = [ABCDE[i] for i in indpr[plr-1]]
    msl = ['u0','v0','w0','expr','expr2']
    matarch[msl]=matarch[msl].astype(str)
    return ABCDE, abc, matarch, plr

# genpolycoeffr(["g","h","alpha"],"z=k0+k1*x+k2*y+k3*x**2+k4*x*y+k5*y**2+k6*x**3+k7*x**2*y+k8*y**2*x+k9*y**3",'alpha',[1,2,4])[2]

def pareq(ste='(x + y*zeta)/(zeta + 1)', **kwargs):
    return eval(ste, kwargs)

prv = ["g","h","alpha"]
rve = [1, 2, 4]
gra = [-44.9573,7.001545,-19.9677]
eqs = ["z=k0+k1*x+k2*y+k3*x**2+k4*x*y+k5*y**2+k6*x**3+k7*x**2*y+k8*y**2*x+k9*y**3","alpha=k0+k1*g+k2*h+k3*Gamma"][0]
sfo = ['g','h','alpha'][2]
test = genpolycoeffr(elem=prv,expr=eqs,solv=sfo,eur=rve)
m1v = math.cos(math.radians(gra[0]))
m2v = math.cos(math.radians(gra[1]))
m3v = math.cos(math.radians(gra[2]))
n1v = math.sin(math.radians(gra[0]))
n2v = math.sin(math.radians(gra[1]))
n3v = math.sin(math.radians(gra[2]))
m1v, m2v, m3v, n1v, n2v, n3v
k0=0.001643394  
k1=1.075717880 
k2=-0.731335532
k3=0.064912642  
k4=0.064007760 
k5=-0.785443920
k6=0.092375583 
k7=-1.271709722  
k8=3.038048882
k9=0.923841730

# I
a1 = pareq(ste=str(test[1][0]),m1=m1v,m2=m2v,m3=m3v,n1=n1v,n2=n2v,n3=n3v)
a2 = pareq(ste=str(test[1][1]),m1=m1v,m2=m2v,m3=m3v,n1=n1v,n2=n2v,n3=n3v)
a3 = pareq(ste=str(test[1][2]),m1=m1v,m2=m2v,m3=m3v,n1=n1v,n2=n2v,n3=n3v)
b1 = pareq(ste=str(test[1][3]),m1=m1v,m2=m2v,m3=m3v,n1=n1v,n2=n2v,n3=n3v)
b2 = pareq(ste=str(test[1][4]),m1=m1v,m2=m2v,m3=m3v,n1=n1v,n2=n2v,n3=n3v)
b3 = pareq(ste=str(test[1][5]),m1=m1v,m2=m2v,m3=m3v,n1=n1v,n2=n2v,n3=n3v)
c1 = pareq(ste=str(test[1][6]),m1=m1v,m2=m2v,m3=m3v,n1=n1v,n2=n2v,n3=n3v)
c2 = pareq(ste=str(test[1][7]),m1=m1v,m2=m2v,m3=m3v,n1=n1v,n2=n2v,n3=n3v)
c3 = pareq(ste=str(test[1][8]),m1=m1v,m2=m2v,m3=m3v,n1=n1v,n2=n2v,n3=n3v)



results = []
for nr in range(20):  # Example: looping from 0 to 19, you can change the range as needed
    lhs = test[2]['d'][nr]
    rhs = test[2]['expr'][nr]
    pareq_result = plsv.pareq(
        str(rhs),
        a1=a1, a2=a2, a3=a3,
        b1=b1, b2=b2, b3=b3,
        c1=c1, c2=c2, c3=c3,
        k0=k0, k1=k1, k2=k2, k3=k3, k4=k4, k5=k5, k6=k6, k7=k7, k8=k8, k9=k9
    )
    # Append the result dictionary directly to the list
    results.append({
        'lhs': lhs,
        'rhs': rhs,
        'pareq_result': pareq_result
    })

df_results = pd.DataFrame(results)

# Output
a1, a2, a3, b1, b2, b3, c1, c2, c3
df_results



