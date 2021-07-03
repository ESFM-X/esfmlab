from sympy import *
import numpy as np
import pandas as pd

class integracion_numerica():
    """
        Aproximación de integrales simples y dobles por los métodos de:
        
        Trapezoidal
            .trapezoidal()
            .trapezoidal_compuesto(particiones)
            .trapecio_compuesto_doble(intervalo2, particiones)
            
        Simpson 1/3
            .simpson1_3()
            .simpson1_3_compuesto(particiones)
            .simpson1_3_compuesto_doble(intervalo2, particiones)
            
        Simpson 3/8
            .simpson3_8()
            .simpson3_8_compuesto(particiones)
            .simpson3_4_compuesto_doble(intervalo2, particiones)
            
        OBS: Cada vez que se ejecute un método nuevo, se tiene que reinstanciar el objeto.
            
        Parámetros
        -----------------------
        limites: list 
            Lista con dos elementos representando los limites de la integral. Si es una
            integral doble, entonces son los limites de la integral de adentro.
        funcion_texto: str
            Representa la función escrita con los operadores de Python.
        
        Atributos
        -----------------------
        a: float
            Limite a de la integral. 
        b: float
            Limite b de la integral. 
        c: float
             Cuando hay dos integrales, representa el límite a de la integral
             de afuera.
        d: float
            Cuando hay dos integrales, representa el límite b de la integral
             de afuera.
             
        solucion: sympy.core.numbers.Float
            Aproximacion a la integral.
        exp: sympy.parse_expr
            Objeto que representa la función simbólica.
        metodo: str
            Nombre del método utilizado para aproximar la solución. Se actualiza
            cada vez que se cambia el método.
        
        aproximado: sympy.core.numbers.Float
            Error aproximado
        estimado: sympy.core.numbers.Float
            Error estimado (es el mismo que cota)
        cota: sympy.core.numbers.Float
            Cota del error (es el mismo que el estimado)
        total: sympy.core.numbers.Float
            Error total
        relativo: sympy.core.numbers.Float
            Error relativo
        verdadero: sympy.core.numbers.Float
            Error verdero
        
        """
    def __init__(self, limites, funcion_texto):
        
        self.a = limites[0]
        self.b = limites[1]
        self.c = None
        self.d = None
       
        self.solucion = None
        self.exp = parse_expr(funcion_texto)

        #self.h = (self.b-self.a)/2
        self.metodo = None
        
        self.aproximado = None
        self.estimado = None
        self.cota = None 
        self.total = None
        self.relativo = None
        self.verdadero = None
        
        self.pasos = []
        
        
    ####----- SIMPLES: ------####  
    def trapezoidal(self):
        x = symbols('x')
        self.solucion =  ((self.b-self.a)/2)*(self.exp.subs(x, self.a) + self.exp.subs(x, self.b))
        self.metodo = "Trapezoidal"

        segunda = integrate(diff(self.exp, x, x), (x, self.a, self.b))/(self.b-self.a)
        self.aproximado = (-((self.b - self.a)**3)/12)*segunda 
        self.estimado =  ((self.b-self.a)**3/12)*self.maximo(3, segunda)
        self.total = (-((self.b - self.a)**3)/12)*diff(self.exp, x, x).subs(x, (self.b-self.a)//2)
        return N(self.solucion)
    
    def simpson1_3(self):
        x = symbols('x')
        h = (self.b-self.a)/2
        self.solucion = h/3 * (self.exp.subs(x,self.a) + 4*self.exp.subs(x,(self.a+self.b)/2) + self.exp.subs(x,self.b) )
        self.metodo = "Simpson 1/3"
        
        tprima = diff(func, x, x, x)
        cuatriprima = diff(func, x, x, x, x)
        
        self.aproximado =  integrate(cuatriprima, (x, self.a,self.b)) 
        self.cota =  ((self.b-self.a)**3/12)*self.maximo(4,tprima)
        
        return N(self.solucion)
    
    def simpson3_8(self):
        h = (self.b-self.a)/3
        x = symbols('x')
        x0 = self.a
        x1 = x0+h
        x2 = x1+h
        x3 = self.b
        
        self.solucion = (x3-x1)*((self.exp.subs(x,x0)+3*self.exp.subs(x,x1)+3*self.exp.subs(x,x2)+self.exp.subs(x,x3)))/8
        self.metodo = "Simpson 3/8"
        
        triprima = diff(fun, x, x, x)
        cuatriprima = diff(fun, x, x, x, x)
        
        self.aproximado = -1*(h**5/90)*((integrate(cuatriprima, (x, x0, x2)))/(x3-x0))
                    
        self.cota = abs((((x3-x0)**5)/6480)*self.maximo(4, triprima))
        
        return self.solucion
    
    ####----- COMPUESTOS: ------####
    
    def trapezoidal_compuesto(self, particiones, errores = True):
        x = symbols('x')
        
        h = (self.b-self.a)/particiones

        self.pasos.append({ 'titulo':'Calcular h', 
                            'procedimiento': '\\( h = \\frac{b-a}{particiones}  \\)',
                            'resultado': '\\( \\Rightarrow  h = \\frac{'+ str(self.b) + ' - ' + str(self.a) +'}{'+str(particiones)+'}  = ' + str(h)+ '\\)'})

        puntos_soporte = np.linspace(self.a + h,self.b - h, particiones-1 )

        self.pasos.append({ 'titulo':'Calcular puntos de soporte', 
                            'procedimiento': 'De h en h desde \\(a\\) hasta \\(b\\)',
                            'resultado': '\\( x_i =  \\)' + ' ' + str(list(puntos_soporte))})

        funcion_evaluada = [self.exp.subs(x,t) for t in puntos_soporte]
        self.pasos.append({ 'titulo':'Evaluar los puntos de soporte en la función', 
                            'procedimiento': '\\( f(x_i) \\)',
                            'resultado': '\\( f(x_i) =  \\)' + ' ' +  str(funcion_evaluada)})

        aprox_inter = sum(funcion_evaluada)
        self.pasos.append({ 'titulo':'Sumar los puntos de soporte evaluados', 
                            'procedimiento': '\\( \\sum f(x_i) \\)',
                            'resultado': '\\( \\sum f(x_i) = \\ \\) ' + str(aprox_inter)})
        fa = self.exp.subs(x,self.a) 
        fb = self.exp.subs(x,self.b)
        self.solucion = h*((1/2)*(fa+ fb) + aprox_inter) 

        self.pasos.append({ 'titulo':'Calcular la aproximación con la fórmula', 
                            'procedimiento': '\\( h\\cdot(\\frac{1}{2} \\cdot (f(a) + f(b)) + \\sum f(x_i) ) \\)',
                            'resultado': f'''\\( \\Rightarrow  {h}\\cdot(\\frac{1}{2} \\cdot (f({self.a}) + f({self.b})) + {aprox_inter} ) \\ \\)
                                            
                                            \\( \\Rightarrow  {h}\\cdot(0.5 \\cdot ({fa} + {fb}) + {aprox_inter} ) = \\ \\)
                                        ''' + str(self.solucion)})

        self.metodo = "Trapezoidal compuesto"
        
        if errores:
            f_biprima = diff(self.exp,x,x)
            
            self.total = abs((-((self.b-self.a)*h**2)/12)*f_biprima.subs(x,(self.b-self.a)/2))
            self.aproximado = (-(h**2)/12)*integrate(f_biprima, (x, self.a, self.b))
            self.cota = (((self.b-self.a)*h**2)/12)*self.maximo(3, f_biprima)
        
        
        return N(self.solucion)
    
    def simpson1_3_compuesto(self, particiones, errores = True):
        x = symbols('x')
        h = (self.b-self.a)/(2*particiones)
        self.pasos.append({ 'titulo':'Calcular h', 
                            'procedimiento': '\\( h = \\frac{b-a}{2 \\cdot particiones}  \\)',
                            'resultado': '\\( \\Rightarrow  h = \\frac{'+ str(self.b) + ' - ' + str(self.a) +'}{2 \\cdot'+str(particiones)+'}  = ' + str(h)+ '\\)'})

        soportes = np.linspace(self.a + h, self.b-h, 2*particiones - 1)
        self.pasos.append({ 'titulo':'Calcular puntos de soporte', 
                            'procedimiento': 'De h en h desde \\(a\\) hasta \\(b\\)',
                            'resultado': '\\( x_i =  \\)' + ' ' + str(list(soportes))})

        cuatri = diff(self.exp, x, x, x, x)

        self.pasos.append({ 'titulo':'Derivar 4 veces \\( \\ f(x) \\)', 
                            'procedimiento': '\\(f^{(i)}(x) = '+latex(diff(self.exp, x)) + ' \\)' + ',  ' +'\\(f^{(ii)}(x) = '+latex(diff(self.exp, x, x)) + ' \\)' + ',  ' +'\\(f^{(iii)}(x) = '+latex(diff(self.exp, x,x,x)) + ' \\)' + ',  ' +'\\(f^{(iv)}(x) = '+latex(cuatri) + ' \\)',
                            'resultado': ''})

        try:
            grado = degree(self.exp, gen = x )
        except PolynomialError:
            Rt = 0
        else:
            if grado > 3:
                Rt = - ((h**5)/90)*cuatri.subs(x,((self.b-self.a)/2))
                self.pasos.append({ 'titulo':'Calcular  \\(\\ R_t \\)', 
                            'procedimiento': 'Debido a que la función es polinómica de grado '+ str(grado) + ', se calcula \\( \\ R_t = \\frac{h^5}{90} \\cdot f^{(iv)} (p) \\), con \\( \\ p \\in (a,b)  \\)',
                            'resultado': '\\( \\Rightarrow  \\ R_t = \\frac{'+str(h)+'^5}{90} \\cdot f^{(iv)} ('+ str((self.b-self.a)/2) +')  \\)' +  '\\( \\  \\Rightarrow  \\ R_t = ' + str((h**5)/90) +'\\cdot' +  str(cuatri.subs(x,((self.b-self.a)/2))) + '  = \\  \\)' + str(Rt)
                            })
            else: 
                Rt = 0

        S1_el = [soportes[i] for i in range(0,2*particiones,2)]
        S2_el = [soportes[i] for i in range(1,2*particiones-1,2)]
        S1_ev = [self.exp.subs(x,i) for i in S1_el]
        S2_ev = [self.exp.subs(x,i) for i in S2_el]
        S1 = sum(S1_ev)
        S2 = sum(S2_ev)

        self.pasos.append({ 'titulo':'Evaluar los puntos de soporte en  \\( \\ f(x) \\)', 
                            'procedimiento': 'Para facilitar cálculos, se divide en dos sumas: \\( \\ S_1 = \\sum_{i=0}^{2\\cdot \\ particiones}x_{2i} \\ \\)  y \\(  \\ S_2 = \\sum_{i=1}^{2\\cdot \\ particiones}x_{2i-1}  \\)',
                            'resultado': 'Puntos de soporte para \\( \\ S_1 \\Rightarrow \\  \\)' + str(S1_el) +'\nEvaluados: '+str(S1_ev) + '. \n\nPara\\( \\ S_2  \\Rightarrow \\ \\)' + str(S2_el) +'.' +' \nEvaluados: '+str(S2_ev)  })
        
        self.pasos.append({ 'titulo':'Calcular  \\( \\ S_1 \\ \\) y \\( \\ S_2 \\)  ', 
                            'procedimiento': 'Como recordatorio  \\( \\ S_1 \\) es la suma de los puntos de soporte en posición par y  \\( \\ S_2 \\) en posición impar. ',
                            'resultado': ' \\( \\ S_1 = \\  \\)' + str(S1) + ' \n\\( \\ S_2 = \\ \\)' + str(S2)  })
        
        fa = self.exp.subs(x, self.a)
        fb = self.exp.subs(x, self.b)

        self.solucion = (h/3)*(fa + 4*S1 + 2*S2 + fb) + Rt

        self.pasos.append({ 'titulo':'Calcular la aproximación con la fórmula', 
                            'procedimiento': '\\( h\\cdot\\frac{1}{3} \\cdot (f(a) + f(b) + 4 \\cdot S1 + 2 \\cdot S2 ) + R_t \\)',
                            'resultado': ' \\( \\Rightarrow \\ ' + str(h) +'\\cdot\\frac{1}{3} \\cdot ('+ str(fa) +' + '+ str(fb) +' + 4 \\cdot '+  str(S1) +'+ 2 \\cdot '+ str(S2) +' ) +'+ str(Rt) +' \\)' +'\n =  ' + str(self.solucion)   })
        
        self.metodo = "Simpson 1/3 compuesto"
        if errores:
            self.total = -(((self.b-self.a)**5)/(180*particiones**4))*cuatri.subs(x,(self.b-self.a)/2) 
            self.aproximado = -((h**4)/180)*(integrate(cuatri,(x, self.a, self.b)))
            self.cota = ((self.b-self.a)*h**4)/180*self.maximo(5,cuatri)
        
        return N(self.solucion)
    
    def simpson3_8_compuesto(self, particiones, errores = False):
        x = symbols('x')
        h = (self.b-self.a)/(3*particiones)
        self.pasos.append({ 'titulo':'Calcular h', 
                            'procedimiento': '\\( h = \\frac{b-a}{2 \\cdot particiones}  \\)',
                            'resultado': '\\( \\Rightarrow  h = \\frac{'+ str(self.b) + ' - ' + str(self.a) +'}{3 \\cdot'+str(particiones)+'}  = ' + str(h)+ '\\)'})

        soportes = np.linspace(self.a+h, self.b-h, 3*particiones - 1)
        self.pasos.append({ 'titulo':'Calcular puntos de soporte', 
                            'procedimiento': 'De h en h desde \\(a\\) hasta \\(b\\)',
                            'resultado': '\\( x_i =  \\)' + ' ' + str(list(soportes))})

        S1_el = [soportes[i] for i in range(0,3*particiones,3)]
        S2_el = [soportes[i] for i in range(1,3*particiones,3)]
        S3_el = [soportes[i] for i in range(2,3*particiones-1,3)]

        S1_ev = [self.exp.subs(x,i) for i in S1_el]
        S2_ev = [self.exp.subs(x,i) for i in S2_el]
        S3_ev = [self.exp.subs(x,i) for i in S3_el]
        
        S1 = sum(S1_ev)# 1,4,7,10
        S2 = sum(S2_ev)# 2,5,8,11
        S3 = sum(S3_ev)# 3,6,9
        #S4 = sum([self.exp.subs(x,soportes[i]) for i in range(3,3*particiones,3)]) 4,

        self.pasos.append({ 'titulo':'Evaluar los puntos de soporte en  \\( \\ f(x) \\)', 
                            'procedimiento': 'Para facilitar cálculos, se divide en dos sumas: \\( \\ S_1 = \\sum_{i=0}^{2\\cdot \\ particiones-1}x_{3i} \\ \\), \\(  \\ S_2 = \\sum_{i=1}^{2\\cdot \\ particiones}x_{3i-1}  \\) y \\(  \\ S_2 = \\sum_{i=2}^{2\\cdot \\ particiones}x_{3i-2}  \\)',
                            'resultado': 'Puntos de soporte para \\( \\ S_1 \\Rightarrow \\  \\)' + str(S1_el) +'\nEvaluados: '+str(S1_ev) + ' \n\nPara\\( \\ S_2  \\Rightarrow \\ \\)' + str(S2_el) +'' +' \nEvaluados: '+str(S2_ev)+ ' \n\nPara\\( \\ S_3  \\Rightarrow \\ \\)' + str(S3_el) +'' +' \nEvaluados: '+str(S3_ev)  })
        
        #print(f'{self.exp.subs(x,self.a)} {self.exp.subs(x,self.b)} 3*{[self.exp.subs(x,soportes[i]) for i in range(0,3*particiones,3)]} 3*{[self.exp.subs(x,soportes[i]) for i in range(1,3*particiones,3)]}  2*{[self.exp.subs(x,soportes[i]) for i in range(2,3*particiones-1,3)]}' )
        self.pasos.append({ 'titulo':'Calcular  \\( \\ S_1 \\ \\),  \\( \\ S_2 \\ \\) y  \\( \\ S_3 \\)  ', 
                            'procedimiento': 'Como recordatorio  \\( \\ S_1 \\) es la suma de los puntos de soporte en posición 0, 3, 4, 7,  10, ...   \\( \\ S_2 \\) en posición 1, 4, 7, 10, ... y  \\( \\ S_3 \\) en posición 2, 5, 8, 11, ... . ',
                            'resultado': ' \\( \\ S_1 = \\  \\)' + str(S1) + ' \n\\( \\ S_2 = \\ \\)' + str(S2) + ' \n\\( \\ S_3 = \\ \\)' + str(S3) })
        
        fa = self.exp.subs(x,self.a)
        fb = self.exp.subs(x,self.b)
        s1_3 = 3*S1
        s2_3 = 3*S2
        s3_2 = 2*S3
        self.solucion = (3*h/8)*(fa + s1_3 + s2_3 + s3_2  + fb)

        self.pasos.append({ 'titulo':'Calcular la aproximación con la fórmula', 
                            'procedimiento': '\\( 3 \\cdot h\\cdot\\frac{1}{8} \\cdot (f(a) + f(b) + 3 \\cdot S1 + 3 \\cdot S2 + 2 \\cdot S3) \\)',
                            'resultado': ' \\( \\Rightarrow \\ 3 \\cdot' + str(h) +'\\cdot\\frac{1}{8} \\cdot ('+ str(fa) +' + '+ str(fb) +' + 3 \\cdot '+  str(S1) +'+ 3\\cdot '+ str(S2) + '\\) \n \\( + 2\\cdot '+ str(S3)  +' \\)' +' =  ' + str(self.solucion)   })
        

        self.metodo = "Simpson 3/8 compuesto"
        if errores:
            cuatriprima = diff(self.exp, x, x, x, x)
            
            self.total = (-((self.b-self.a)/80)*h**4)*cuatriprima.subs(x,(self.b-self.a)/2)
            self.aproximado = (-(h**4/80))*(integrate(cuatriprima, (x, self.a, self.b)))
            self.cota = ((self.b-self.a)*h**4)/80*self.maximo(5,cuatriprima)
        
        return  self.solucion
    
    ####----- DOBLES: ------####
    def trapecio_compuesto_doble(self, intervalo2, particiones):
        """
        Recibe los limites de la integral de afuera.
        """
        x = symbols('x')
        y = symbols('y')
        try: 
            a = float(self.a)
            b = float(self.b)
        except:
            bandera = 0
        else:
            bandera = 1
        
        if bandera:
            h = (self.b-self.a)/particiones
            self.pasos.append({ 'titulo':'Calcular \\( \\ h_x \\ \\)', 
                            'procedimiento': '\\( h_x = \\frac{b-a}{ particiones}  \\)',
                            'resultado': '\\( \\Rightarrow  h_x = \\frac{'+ str(self.b) + ' - ' + str(self.a) +'}{'+str(particiones)+'}  = ' + str(h)+ '\\)'})


            puntos_soporte = np.linspace(self.a + h,self.b - h, particiones-1)

            self.pasos.append({ 'titulo':'Calcular puntos de soporte', 
                            'procedimiento': 'De \\( \\ h_x \\ \\) en \\( \\ h_x \\ \\) desde \\(a\\) hasta \\(b\\)',
                            'resultado': '\\( x_i =  \\)' + ' ' + str(list(puntos_soporte))})

            funcion_evaluada = [self.exp.subs(x,t) for t in puntos_soporte]
            funcion_evaluada_latex = [' \\(  ' + latex(t) + ' \\) ' for t in funcion_evaluada]
            self.pasos.append({ 'titulo':'Evaluar los puntos de soporte en  \\( \\ f(x,y) \\)', 
                            'procedimiento': '\\( f(x_i, y) \\)',
                            'resultado': '\\( f(x_i, y) = \\ \\)' + '[ '+ (', '.join(funcion_evaluada_latex)) + ']' })
        
            aprox_inter = sum(funcion_evaluada)
            self.pasos.append({ 'titulo':'Calcular la suma de los puntos de soporte evaluados', 
                            'procedimiento': '\\( \\sum f(x_i, y) \\)',
                            'resultado': ' = \\( ' + latex(aprox_inter) + '\\)' })

            fa = self.exp.subs(x,self.a)
            fb = self.exp.subs(x,self.b)
    
            primera_integral =  h*((1/2)*(fa + fb) + aprox_inter) 

            self.pasos.append({ 'titulo':'Calcular \\( \\ g(y) \\ \\)con la fórmula', 
                            'procedimiento': '\\( h_x\\cdot(\\frac{1}{2} \\cdot (f(a,y) + f(b,y)) + \\sum f(x_i, y)) \\)',
                            'resultado': '\\( '+ str(h) + ' \\cdot(\\frac{1}{2} \\cdot ('+latex(fa) +' + ' + latex(fb) +') + ' + latex(aprox_inter)  +') = ' + latex(simplify(primera_integral)) + ' \\)' + '\n \\( \\therefore g(y) = '+  latex(simplify(primera_integral)) + ' \\)'   })
        

            c = intervalo2[0]
            d = intervalo2[1]
            h = (d-c)/particiones
            self.pasos.append({ 'titulo':'Calcular \\( \\ h_y \\ \\)', 
                            'procedimiento': '\\( h_y = \\frac{d-c}{ particiones}  \\)',
                            'resultado': '\\( \\Rightarrow  h_y = \\frac{'+ str(d) + ' - ' + str(c) +'}{'+str(particiones)+'}  = ' + str(h)+ '\\)'})


            self.c = c
            self.d = d
            puntos_soporte = np.linspace(c + h, d - h, particiones-1)

            self.pasos.append({ 'titulo':'Calcular puntos de soporte', 
                            'procedimiento': 'De \\( \\ h_y \\ \\) en \\( \\ h_y \\ \\) desde \\(c\\) hasta \\(d\\)',
                            'resultado': '\\( y_i =  \\)' + ' ' + str(list(puntos_soporte))})

            funcion_evaluada = [primera_integral.subs(y,t) for t in puntos_soporte]
            funcion_evaluada_latex =  [' \\(  ' + latex(t) + ' \\) ' for t in funcion_evaluada]
            self.pasos.append({ 'titulo':'Evaluar los puntos de soporte en  \\( \\ g(y) \\)', 
                            'procedimiento': '\\( g(y_i) \\)',
                            'resultado': '\\( g(y_i) = \\ \\)' + '[ '+ (', '.join(funcion_evaluada_latex)) + ']' })
        
            aprox_inter = sum(funcion_evaluada)

            self.pasos.append({ 'titulo':'Calcular la suma de los puntos de soporte evaluados', 
                            'procedimiento': '\\( \\sum g(y_i) \\)',
                            'resultado': ' = \\( ' + latex(aprox_inter) + '\\)' })

            gc = primera_integral.subs(y,c)
            gd = primera_integral.subs(y,d)
            self.solucion =  h*((1/2)*(gc + gd) + aprox_inter) 

            self.pasos.append({ 'titulo':'Aproximar la integral con la fórmula', 
                            'procedimiento': '\\( h_y\\cdot(\\frac{1}{2} \\cdot (g(c) + g(d)) + \\sum g(y_i)) \\)',
                            'resultado': '\\( '+ str(h) + ' \\cdot(\\frac{1}{2} \\cdot ('+ str(gc) +' + ' + str(gd) +') + ' + str(aprox_inter)  +') = ' + str(self.solucion) +' \\)'   })
        
            self.metodo = "Trapezoidal compuesto doble numérico"

            f_biprima = diff(self.exp,x,x)

            self.total = abs((-((self.b-self.a)*h**2)/12)*f_biprima.subs(x,(self.b-self.a)/2))
            self.aproximado = (-(h**2)/12)*integrate(f_biprima, (x, self.a, self.b))
            self.cota = (((self.b-self.a)*h**2)/12)*self.maximo(3, f_biprima)

            return N(self.solucion)
        else:
            c = intervalo2[0]
            d = intervalo2[1]
            h = (d-c)/particiones
            self.pasos.append({ 'titulo':'Calcular \\( \\ h_x \\ \\)', 
                            'procedimiento': '\\( h_x = \\frac{d-c}{ particiones}  \\)',
                            'resultado': '\\( \\Rightarrow  h_x = \\frac{'+ str(d) + ' - ' + str(c) +'}{'+str(particiones)+'}  = ' + str(h)+ '\\)'})
            a = parse_expr(self.a)
            b = parse_expr(self.b)
            self.a = a
            self.b = b
            self.c = c
            self.d = d
            puntos_soporte = np.linspace(c + h,d - h, particiones-1)

            self.pasos.append({ 'titulo':'Calcular puntos de soporte', 
                            'procedimiento': 'De \\( \\ h_x \\ \\) en \\( \\ h_x \\ \\) desde \\(c\\) hasta \\(d\\)',
                            'resultado': '\\( x_i =  \\)' + ' ' + str(list(puntos_soporte))})
            funcion_evaluada = []
            procedimiento = []
            for t in puntos_soporte:
                aa = float(a.subs(x, t))
                bb = float(b.subs(x, t))
                expresion =  self.exp.subs(x,t)

                integral = integracion_numerica([aa, bb], str(expresion.subs(y, x)) )
                aproximacion = integral.trapezoidal_compuesto(particiones, errores = False)

                procedimiento.append({'titulo':'Calcular evaluando el punto de soporte \\( \\ '+ str(t) +' \\ \\) \\( \\ \\Rightarrow \\int_{'+ str(aa)+ '}^{'+ str(bb)+ '} '+ latex(expresion) +'\\ dy\\ \\) con el método Trapezoidal de '+ str(particiones) + ' particiones', 
                                      'procedimiento': integral.pasos,
                                      'resultado': str(aproximacion) 
                                    })
                funcion_evaluada.append(aproximacion)
            self.pasos.append({ 'titulo':'Calcular las integrales evaluando los puntos de soporte en los límites de la integral y en la función', 
                            'procedimiento2': procedimiento,
                            'resultado': '\\( G(x_i) =  ' + str(funcion_evaluada) + '\\)' })

            #funcion_evaluada = [integracion_numerica([float(a.subs(x, t)), float(b.subs(x, t))], str(self.exp.subs(x,t).subs(y, x))).trapezoidal_compuesto(particiones) for t in puntos_soporte]
            aprox_inter = sum(funcion_evaluada)

            self.pasos.append({ 'titulo':'Calcular la suma de las integrales evaluadas', 
                            'procedimiento': '\\( \\sum G(x_i) \\)',
                            'resultado': '\\( \\sum G(x_i) =  ' + str(aprox_inter) + '\\)' })

            aa1 = float(a.subs(x, c))
            bb1 = float(b.subs(x,c))
            expression1 = self.exp.subs(x,c)

            aa2 = float(a.subs(x, d))
            bb2 = float(b.subs(x, d))
            expression2 = self.exp.subs(x,d)

            procedimiento2 = []
            integral_c = integracion_numerica([aa1, bb1], str(expression1.subs(y, x)))
            aproximacionc = integral_c.trapezoidal_compuesto(particiones, errores = False) 
            procedimiento2.append({'titulo':'Calcular evaluando el punto  \\( \\ '+ str(c) +' \\ \\) \\( \\ \\Rightarrow \\int_{'+ str(aa1)+ '}^{'+ str(bb1)+ '} '+ latex(expression1) +'\\ dy\\ \\) con el método Trapezoidal de '+ str(particiones) + ' particiones', 
                                      'procedimiento': integral_c.pasos,
                                      'resultado': str(aproximacionc) 
                                    })

            integral_d = integracion_numerica([aa2, bb2], str(expression2.subs(y, x)))
            aproximaciond = integral_d.trapezoidal_compuesto(particiones, errores = False)
            procedimiento2.append({'titulo':'Calcular evaluando el punto  \\( \\ '+ str(d) +' \\ \\) \\( \\ \\Rightarrow \\int_{'+ str(aa2)+ '}^{'+ str(bb2)+ '} '+ latex(expression2) +'\\ dy\\ \\) con el método Trapezoidal de '+ str(particiones) + ' particiones', 
                                      'procedimiento': integral_d.pasos,
                                      'resultado': str(aproximaciond) 
                                    })

            self.pasos.append({ 'titulo':'Calcular las integrales evaluando los límites y la función con \\( \\ c \\ \\) y \\( \\ d \\ \\)', 
                            'procedimiento2': procedimiento2,
                            'resultado': '\\( G(c) =  ' + str(aproximacionc) + '\\)' + '\\( \\ \\ G(d) =  ' + str(aproximaciond) + '\\)'  })

            self.solucion =  h*((1/2)*(aproximacionc + aproximaciond)  + aprox_inter) 
            
            self.pasos.append({ 'titulo':'Aproximar la integral con la fórmula', 
                            'procedimiento': '\\( h_x\\cdot(\\frac{1}{2} \\cdot(G(c) + G(d) ) + \sum G(x_i)  \\) ',
                            'resultado': '\\( \\Rightarrow '+ str(h)+'\\cdot(\\frac{1}{2} \\cdot ('+ str(aproximacionc) +'  + '+ str(aproximaciond) +') +'+ str(aprox_inter) +'  = \\ ' + str(self.solucion)+' \\)'  })

            self.metodo = "Trapezoidal compuesto doble"
            return N(self.solucion)
        
    def simpson1_3_compuesto_doble(self, intervalo2, particiones):
        """
        Recibe los limites de la integral de afuera.
        """
        x = symbols('x')
        y = symbols('y')

        try: 
            a = float(self.a)
            b = float(self.b)
        except:
            bandera = 0
        else:
            bandera = 1
        
        if bandera:
            h = (self.b-self.a)/(2*particiones)
            self.pasos.append({ 'titulo':'Calcular \\( \\ h_x \\ \\)', 
                            'procedimiento': '\\( h_x = \\frac{b-a}{ 2 \\cdot particiones}  \\)',
                            'resultado': '\\( \\Rightarrow  h_x = \\frac{'+ str(self.b) + ' - ' + str(self.a) +'}{2 \\cdot'+str(particiones)+'}  = ' + str(h)+ '\\)'})
            soportes = np.linspace(self.a + h, self.b-h, 2*particiones - 1)
            self.pasos.append({ 'titulo':'Calcular puntos de soporte', 
                            'procedimiento': 'De \\( \\ h_x \\ \\) en \\( \\ h_x \\ \\) desde \\(a\\) hasta \\(b\\)',
                            'resultado': '\\( x_i =  \\)' + ' ' + str(list(soportes))})

            S1_el = [soportes[i] for i in range(0,2*particiones,2)]
            S2_el = [soportes[i] for i in range(1,2*particiones-1,2)]
            S1_ev = [self.exp.subs(x,i) for i in S1_el]
            S2_ev = [self.exp.subs(x,i) for i in S2_el]
            S1 = sum(S1_ev)
            S2 = sum(S2_ev)
            
            self.pasos.append({ 'titulo':'Evaluar los puntos de soporte en  \\( \\ f(x,y) \\)', 
                            'procedimiento': 'Para facilitar cálculos, se divide en dos sumas: \\( \\ S_1 = \\sum_{i=0}^{2\\cdot \\ particiones}f(x_{2i},y) \\ \\)  y \\(  \\ S_2 = \\sum_{i=1}^{2\\cdot \\ particiones}f(x_{2i-1},y)  \\)',
                            'resultado': 'Puntos de soporte para \\( \\ S_1 \\Rightarrow \\  ' + latex(S1_el) +'\\)\nEvaluados: \\( '+latex(S1_ev) + '\\). \n\nPara\\( \\ S_2  \\Rightarrow \\ ' + latex(S2_el) +'\\)' +' \nEvaluados: \\('+latex(S2_ev)+'\\)'  })
            
            self.pasos.append({ 'titulo':'Calcular \\( \\ S_1 \\ \\) y \\( \\ S_2 \\)  ', 
                            'procedimiento': 'Como recordatorio  \\( \\ S_1 \\) es la suma de los puntos de soporte evaluados en posición par y  \\( \\ S_2 \\) en posición impar. ',
                            'resultado': ' \\( \\ S_1 = \\ ' + latex(S1) + ' \\) \n\\( \\ S_2 = \\ ' + latex(S2) +'\\)' })

            fa  = self.exp.subs(x, self.a)
            fb = self.exp.subs(x, self.b)
            primera_integral =  (h/3)*(fa + 4*S1 + 2*S2 + fb) #+Rt
            self.pasos.append({ 'titulo':'Calcular \\( \\ g(y) \\ \\)con la fórmula', 
                            'procedimiento': '\\( h\\cdot\\frac{1}{3} \\cdot (f(a,y) + f(b,y) + 4 \\cdot S1 + 2 \\cdot S2 )  \\)',
                            'resultado': ' \\( \\Rightarrow \\ ' + str(h) +'\\cdot\\frac{1}{3} \\cdot ('+ latex(fa) +' + '+ latex(fb) +' + 4 \\cdot '+  latex(S1) +'+ 2 \\cdot '+ latex(S2) +' )  \\)' +'\n \\(= \\ ' + latex(primera_integral)+'\\)'    })
            c = intervalo2[0]
            d = intervalo2[1]
            h = (d-c)/(2*particiones)
            self.pasos.append({ 'titulo':'Calcular \\( \\ h_y \\ \\)', 
                            'procedimiento': '\\( h_y = \\frac{d-c}{ 2 \\cdot particiones}  \\)',
                            'resultado': '\\( \\Rightarrow  h_y = \\frac{'+ str(d) + ' - ' + str(c) +'}{2 \\cdot'+str(particiones)+'}  = ' + str(h)+ '\\)'})
            

            self.c = c
            self.d = d
            soportes = np.linspace(c + h, d - h, 2*particiones - 1)
            self.pasos.append({ 'titulo':'Calcular puntos de soporte', 
                            'procedimiento': 'De \\( \\ h_y \\ \\) en \\( \\ h_y \\ \\) desde \\(c\\) hasta \\(d\\)',
                            'resultado': '\\( y_i =  \\)' + ' ' + str(list(soportes))})

            S1_el = [soportes[i] for i in range(0,2*particiones,2)]
            S2_el = [soportes[i] for i in range(1,2*particiones-1,2)]
            S1_ev = [primera_integral.subs(y,i) for i in S1_el]
            S2_ev = [primera_integral.subs(y,i) for i in S2_el]
            S1 = sum(S1_ev)
            S2 = sum(S2_ev)

            self.pasos.append({ 'titulo':'Evaluar los puntos de soporte en  \\( \\ g(y) \\)', 
                            'procedimiento': 'Para facilitar cálculos, se divide en dos sumas: \\( \\ S_1 = \\sum_{i=0}^{2\\cdot \\ particiones}g(y_{2i}) \\ \\)'+  '  y ' +'\\(  \\ S_2 = \\sum_{i=1}^{2\\cdot \\ particiones}g(y_{2i-1})  \\)',
                            'resultado': 'Puntos de soporte para \\( \\ S_1 \\Rightarrow \\  ' + latex(S1_el) +'\\)\nEvaluados: \\( '+latex(S1_ev) + '\\). \n\nPara\\( \\ S_2  \\Rightarrow \\ ' + latex(S2_el) +'\\)' +' \nEvaluados: \\('+latex(S2_ev)+'\\)'  })
            
            self.pasos.append({ 'titulo':'Calcular  \\( \\ S_1 \\ \\) y \\( \\ S_2 \\)  ', 
                            'procedimiento': 'Como recordatorio  \\( \\ S_1 \\) es la suma de los puntos de soporte evaluados en posición par y  \\( \\ S_2 \\) en posición impar. ',
                            'resultado': ' \\( \\ S_1 = \\ ' + latex(S1) + ' \\) \n\\( \\ S_2 = \\ ' + latex(S2) +'\\)' })

            gc = primera_integral.subs(y, c)
            gd = primera_integral.subs(y, d)           
            self.solucion =  (h/3)*(gc + 4*S1 + 2*S2 + gd) #+Rt
            self.pasos.append({ 'titulo':'Aproximar la integral con la fórmula', 
                            'procedimiento': '\\( h\\cdot\\frac{1}{3} \\cdot (g(c) + g(d) + 4 \\cdot S1 + 2 \\cdot S2 )  \\)',
                            'resultado': ' \\( \\Rightarrow \\ ' + str(h) +'\\cdot\\frac{1}{3} \\cdot ('+ latex(gc) +' + '+ latex(gd) +' + 4 \\cdot '+  latex(S1) +'+ 2 \\cdot '+ latex(S2) +' )  \\)' +'\n \\(= \\ ' + latex(self.solucion)+'\\)'    })
            
            self.metodo = "Simpson 1/3 compuesto doble numérico"

            return N(self.solucion)
            
        else:
            c = intervalo2[0]
            d = intervalo2[1]
            h = (d-c)/(2*particiones)

            self.pasos.append({ 'titulo':'Calcular \\( \\ h_x \\ \\)', 
                            'procedimiento': '\\( h_x = \\frac{d-c}{2\\cdot particiones}  \\)',
                            'resultado': '\\( \\Rightarrow  h_x = \\frac{'+ str(d) + ' - ' + str(c) +'}{ 2 \\cdot '+str(particiones)+'}  = ' + str(h)+ '\\)'})


            a = parse_expr(self.a)
            b = parse_expr(self.b)
            
            self.a = a
            self.b = b
            self.c = c
            self.d = d
            
            soportes = np.linspace(c + h, d - h, 2*particiones - 1)
            
            self.pasos.append({ 'titulo':'Calcular puntos de soporte', 
                            'procedimiento': 'De \\( \\ h_x \\ \\) en \\( \\ h_x \\ \\) desde \\(c\\) hasta \\(d\\)',
                            'resultado': '\\( x_i =  \\)' + ' ' + str(list(soportes))})
            
            ##### s1
            funcion_evaluada = []
            procedimiento = []
            puntos_soporteS1 = [soportes[t] for t in range(0,2*particiones,2)]
            for t in puntos_soporteS1:
                aa = float(a.subs(x, t))
                bb = float(b.subs(x, t))
                expresion =  self.exp.subs(x,t)

                integral = integracion_numerica([aa, bb], str(expresion.subs(y, x)) )
                aproximacion = integral.simpson1_3_compuesto(particiones, errores = False)

                procedimiento.append({'titulo':'Calcular evaluando el punto de soporte \\( \\ '+ str(t) +' \\ \\) \\( \\ \\Rightarrow \\int_{'+ str(aa)+ '}^{'+ str(bb)+ '} '+ latex(expresion) +'\\ dy\\ \\) con el método Simpson 1/3 de '+ str(particiones) + ' particiones', 
                                      'procedimiento': integral.pasos,
                                      'resultado': str(aproximacion) 
                                    })
                funcion_evaluada.append(aproximacion)
            S1 = sum(funcion_evaluada)

            self.pasos.append({ 'titulo':'Calcular S1 con las integrales evaluando los puntos de soporte pares en los límites de la integral y en la función', 
                            'procedimiento2': procedimiento,
                            'resultado': '\\( S1 = \\sum \\ ' + str(funcion_evaluada) + ' \\ = \\ '+ str(S1) + ' \\)' })
            ##### s2
            funcion_evaluada = []
            procedimiento = []
            puntos_soporteS1 = [soportes[t] for t in range(1,2*particiones-1,2)]
            for t in puntos_soporteS1:
                aa = float(a.subs(x, t))
                bb = float(b.subs(x, t))
                expresion =  self.exp.subs(x,t)

                integral = integracion_numerica([aa, bb], str(expresion.subs(y, x)) )
                aproximacion = integral.simpson1_3_compuesto(particiones, errores = False)

                procedimiento.append({'titulo':'Calcular evaluando el punto de soporte \\( \\ '+ str(t) +' \\ \\) \\( \\ \\Rightarrow \\int_{'+ str(aa)+ '}^{'+ str(bb)+ '} '+ latex(expresion) +'\\ dy\\ \\) con el método Simpson 1/3 de '+ str(particiones) + ' particiones', 
                                      'procedimiento': integral.pasos,
                                      'resultado': str(aproximacion) 
                                    })
                funcion_evaluada.append(aproximacion)
            S2 = sum(funcion_evaluada)
            
            self.pasos.append({ 'titulo':'Calcular S2 con las integrales evaluando los puntos de soporte impares en los límites de la integral y en la función', 
                            'procedimiento2': procedimiento,
                            'resultado': '\\( S2 = \\sum \\ ' + str(funcion_evaluada) + ' \\ = \\ '+ str(S2) + ' \\)' })
            
            aa1 = float(a.subs(x, c))
            bb1 = float(b.subs(x,c))
            expression1 = self.exp.subs(x,c)

            aa2 = float(a.subs(x, d))
            bb2 = float(b.subs(x, d))
            expression2 = self.exp.subs(x,d)

            procedimiento2 = []
            integral_gc = integracion_numerica([aa1, bb1], str(expression1.subs(y, x)))
            gc = integral_gc.simpson1_3_compuesto(particiones, errores = False)
            procedimiento2.append({'titulo':'Calcular evaluando el punto  \\( \\ '+ str(c) +' \\ \\) \\( \\ \\Rightarrow \\int_{'+ str(aa1)+ '}^{'+ str(bb1)+ '} '+ latex(expression1) +'\\ dy\\ \\) con el método Simpson 1/3 de '+ str(particiones) + ' particiones', 
                                      'procedimiento': integral_gc.pasos,
                                      'resultado': str(gc) 
                                    })
            integral_gd = integracion_numerica([aa2, bb2], str(expression2.subs(y, x)))
            gd = integral_gd.simpson1_3_compuesto(particiones, errores = False)
            procedimiento2.append({'titulo':'Calcular evaluando el punto  \\( \\ '+ str(d) +' \\ \\) \\( \\ \\Rightarrow \\int_{'+ str(aa2)+ '}^{'+ str(bb1)+ '} '+ latex(expression2) +'\\ dy\\ \\) con el método Simpson 1/3 de '+ str(particiones) + ' particiones', 
                                      'procedimiento': integral_gd.pasos,
                                      'resultado': str(gd) 
                                    })

            self.pasos.append({ 'titulo':'Calcular las integrales evaluando los límites y la función con \\( \\ c \\ \\) y \\( \\ d \\ \\) ', 
                            'procedimiento2':procedimiento2 ,
                            'resultado': f'\\(G(c) = {gc} \\ \\ G(d) = {gd} \\)'  })


            self.solucion =  (h/3)*(gc + 4*S1 + 2*S2 + gd ) #+Rt
            
            self.pasos.append({ 'titulo':'Calcular la aproximación con la fórmula', 
                            'procedimiento': '\\( h\\cdot\\frac{1}{3} \\cdot (G(c) + G(d) + 4 \\cdot S1 + 2 \\cdot S2 )  \\)',
                            'resultado': ' \\( \\Rightarrow \\ ' + str(h) +'\\cdot\\frac{1}{3} \\cdot ('+ str(gc) +' + '+ str(gd) +' + 4 \\cdot '+  str(S1) +'+ 2 \\cdot '+ str(S2) +' )  \\)' +'\n =  ' + str(self.solucion)   })

            self.metodo = "Simpson 1/3 compuesto doble"
            return N(self.solucion)            
        
    def simpson3_8_compuesto_doble(self, intervalo2, particiones):
        """
        Recibe los limites de la integral de afuera.
        """
        x = symbols('x')
        y = symbols('y')
        try: 
            a = float(self.a)
            b = float(self.b)
        except:
            bandera = 0
        else:
            bandera = 1
        
        if bandera:
           
            h = (self.b-self.a)/(3*particiones)
            self.pasos.append({ 'titulo':'Calcular \\( \\ h_x \\ \\)', 
                            'procedimiento': '\\( h_x = \\frac{b-a}{ 3 \\cdot particiones}  \\)',
                            'resultado': '\\( \\Rightarrow  h_x = \\frac{'+ str(self.b) + ' - ' + str(self.a) +'}{3 \\cdot'+str(particiones)+'}  = ' + str(h)+ '\\)'})
            soportes = np.linspace(self.a+h, self.b-h, 3*particiones - 1)
            self.pasos.append({ 'titulo':'Calcular puntos de soporte', 
                            'procedimiento': 'De \\( \\ h_x \\ \\) en \\( \\ h_x \\ \\) desde \\(a\\) hasta \\(b\\)',
                            'resultado': '\\( x_i =  \\)' + ' ' + str(list(soportes))})

            S1_el = [soportes[i] for i in range(0,3*particiones,3)]
            S2_el = [soportes[i] for i in range(1,3*particiones,3)]
            S3_el = [soportes[i] for i in range(2,3*particiones-1,3)]
            S1_ev = [self.exp.subs(x,i) for i in S1_el]
            S2_ev = [self.exp.subs(x,i) for i in S2_el]
            S3_ev = [self.exp.subs(x,i) for i in S3_el]
            S1 = sum(S1_ev)
            S2 = sum(S2_ev)
            S3 = sum(S3_ev)

            self.pasos.append({ 'titulo':'Evaluar los puntos de soporte en  \\( \\ f(x,y) \\)', 
                            'procedimiento': 'Para facilitar cálculos, se divide en tres sumas: \\( \\ S_1 = \\sum_{i=0}^{2\\cdot \\ particiones}f(x_{3i},y) \\ \\), \\(  \\ S_2 = \\sum_{i=1}^{2\\cdot \\ particiones}f(x_{3i-1},y)  \\)  y \\(  \\ S_3 = \\sum_{i=2}^{2\\cdot \\ particiones}f(x_{3i-3},y)  \\)',
                            'resultado': 'Puntos de soporte para \\( \\ S_1 \\Rightarrow \\  ' + latex(S1_el) +'\\)\nEvaluados: \\( '+latex(S1_ev) + '\\). \n\nPara\\( \\ S_2  \\Rightarrow \\ ' + latex(S2_el) +'\\)' +' \nEvaluados: \\('+latex(S2_ev)+'\\) \n\nPara\\( \\ S_3  \\Rightarrow \\ ' + latex(S3_el) +'\\)' +' \nEvaluados: \\('+latex(S3_ev)+'\\)' })
            
            self.pasos.append({ 'titulo':'Calcular \\( \\ S_1 \\ \\), \\( \\ S_2 \\  \\) y \\( \\ S_3 \\)   ', 
                            'procedimiento': 'Como recordatorio  \\( \\ S_1 \\) es la suma de los puntos de soporte con posiciones de 3 en 3 empezando en 0,  \\( \\ S_2 \\) empezando en 1 y \\( \\ S_2 \\) empezando en 2. ',
                            'resultado': ' \\( \\ S_1 = \\ ' + latex(S1) + ' \\) \n\\( \\ S_2 = \\ ' + latex(S2) +'\\) \n \\( \\ S_3 = \\ ' + latex(S3) +'\\)' })

            fa  = self.exp.subs(x, self.a)
            fb = self.exp.subs(x, self.b)
            primera_integral =  (3*h/8)*(fa + 3*S1 + 3*S2 + 2*S3  + fb)
            self.pasos.append({ 'titulo':'Calcular \\( \\ g(y) \\ \\)con la fórmula', 
                            'procedimiento': '\\( h\\cdot\\frac{3}{8} \\cdot (f(a,y) + f(b,y) + 3 \\cdot S1 + 3 \\cdot S2 + 2 \\cdot S2 )  \\)',
                            'resultado': ' \\( \\Rightarrow \\ ' + str(h) +'\\cdot\\frac{1}{3} \\cdot ('+ latex(fa) +' + '+ latex(fb) +' + 4 \\cdot '+  latex(S1) +'+ 2 \\cdot '+ latex(S2) +' )  \\)' +'\n \\(= \\ ' + latex(primera_integral)+'\\)'    })
            c = intervalo2[0]
            d = intervalo2[1]
            h = (d-c)/(3*particiones)
            self.pasos.append({ 'titulo':'Calcular \\( \\ h_y \\ \\)', 
                            'procedimiento': '\\( h_y = \\frac{d-c}{ 3 \\cdot particiones}  \\)',
                            'resultado': '\\( \\Rightarrow  h_y = \\frac{'+ str(d) + ' - ' + str(c) +'}{3 \\cdot'+str(particiones)+'}  = ' + str(h)+ '\\)'})
            self.c = c
            self.d = d
            soportes = np.linspace(c + h, d - h, 3*particiones - 1)
            self.pasos.append({ 'titulo':'Calcular puntos de soporte', 
                            'procedimiento': 'De \\( \\ h_y \\ \\) en \\( \\ h_y \\ \\) desde \\(c\\) hasta \\(d\\)',
                            'resultado': '\\( y_i =  \\)' + ' ' + str(list(soportes))})

            S1_el = [soportes[i] for i in range(0,3*particiones,3)]
            S2_el = [soportes[i] for i in range(1,3*particiones-1,3)]############################################################# -1?
            S3_el = [soportes[i] for i in range(2,3*particiones-1,3)]
            S1_ev = [primera_integral.subs(y,i) for i in S1_el]
            S2_ev = [primera_integral.subs(y,i) for i in S2_el]
            S3_ev = [primera_integral.subs(y,i) for i in S3_el]
            S1 = sum(S1_ev)
            S2 = sum(S2_ev)
            S3 = sum(S3_ev)

            self.pasos.append({ 'titulo':'Evaluar los puntos de soporte en  \\( \\ g(y) \\)', 
                            'procedimiento': 'Para facilitar cálculos, se divide en tres sumas: \\( \\ S_1 = \\sum_{i=0}^{2\\cdot \\ particiones}g(y_{3i}) \\ \\), \\(  \\ S_2 = \\sum_{i=1}^{2\\cdot \\ particiones}g(y_{3i-1})  \\)  y \\(  \\ S_3 = \\sum_{i=2}^{2\\cdot \\ particiones}g(y_{3i-3})  \\)',
                            'resultado': 'Puntos de soporte para \\( \\ S_1 \\Rightarrow \\  ' + latex(S1_el) +'\\)\nEvaluados: \\( '+latex(S1_ev) + '\\). \n\nPara\\( \\ S_2  \\Rightarrow \\ ' + latex(S2_el) +'\\)' +' \nEvaluados: \\('+latex(S2_ev)+'\\) \n\nPara\\( \\ S_3  \\Rightarrow \\ ' + latex(S3_el) +'\\)' +' \nEvaluados: \\('+latex(S3_ev)+'\\)' })
            
            self.pasos.append({ 'titulo':'Calcular \\( \\ S_1 \\ \\), \\( \\ S_2 \\  \\) y \\( \\ S_3 \\)   ', 
                            'procedimiento': 'Como recordatorio  \\( \\ S_1 \\) es la suma de los puntos de soporte evaluados con posiciones de 3 en 3 empezando en 0,  \\( \\ S_2 \\) empezando en 1 y \\( \\ S_2 \\) empezando en 2. ',
                            'resultado': ' \\( \\ S_1 = \\ ' + latex(S1) + ' \\) \n\\( \\ S_2 = \\ ' + latex(S2) +'\\) \n \\( \\ S_3 = \\ ' + latex(S3) +'\\)' })
            fa = primera_integral.subs(y,c)
            fb = primera_integral.subs(y,d)
            self.solucion =  (3*h/8)*(fa + 3*S1 + 3*S2 + 2*S3  + fb)
            self.pasos.append({ 'titulo':'Calcular la aproximación con la fórmula', 
                            'procedimiento': '\\( 3 \\cdot h\\cdot\\frac{1}{8} \\cdot (g(c) + g(d) + 3 \\cdot S1 + 3 \\cdot S2 + 2 \\cdot S3) \\)',
                            'resultado': ' \\( \\Rightarrow \\ 3 \\cdot' + str(h) +'\\cdot\\frac{1}{8} \\cdot ('+ str(fa) +' + '+ str(fb) +' + 3 \\cdot '+  str(S1) +'+ 3\\cdot '+ str(S2) + '\\) \n \\( + 2\\cdot '+ str(S3)  +' =  ' + str(self.solucion)  +' \\)'   })
            self.metodo = "Simpson 3/8 compuesto doble numérico"

            return N(self.solucion)
            
        else:
            c = intervalo2[0]
            d = intervalo2[1]
            h = (d-c)/(3*particiones)
            self.pasos.append({ 'titulo':'Calcular \\( \\ h_x \\ \\)', 
                            'procedimiento': '\\( h_x = \\frac{d-c}{3\\cdot particiones}  \\)',
                            'resultado': '\\( \\Rightarrow  h_x = \\frac{'+ str(d) + ' - ' + str(c) +'}{ 3 \\cdot '+str(particiones)+'}  = ' + str(h)+ '\\)'})

            a = parse_expr(self.a)
            b = parse_expr(self.b)
            
            self.a = a
            self.b = b
            self.c = c
            self.d = d
            
            soportes = np.linspace(c + h, d - h, 3*particiones - 1)
            self.pasos.append({ 'titulo':'Calcular puntos de soporte', 
                            'procedimiento': 'De \\( \\ h_x \\ \\) en \\( \\ h_x \\ \\) desde \\(c\\) hasta \\(d\\)',
                            'resultado': '\\( x_i =  \\)' + ' ' + str(list(soportes))})

            ##### s1
            funcion_evaluada = []
            procedimiento = []
            puntos_soporteS1 = [soportes[t] for t in range(0,3*particiones,3)]
            for t in puntos_soporteS1:
                aa = float(a.subs(x, t))
                bb = float(b.subs(x, t))
                expresion =  self.exp.subs(x,t)

                integral = integracion_numerica([aa, bb], str(expresion.subs(y, x)) )
                aproximacion = integral.simpson3_8_compuesto(particiones, errores = False)

                procedimiento.append({'titulo':'Calcular evaluando el punto de soporte \\( \\ '+ str(t) +' \\ \\) \\( \\ \\Rightarrow \\int_{'+ str(aa)+ '}^{'+ str(bb)+ '} '+ latex(expresion) +'\\ dy\\ \\) con el método Simpson 3/8 de '+ str(particiones) + ' particiones', 
                                      'procedimiento': integral.pasos,
                                      'resultado': str(aproximacion) 
                                    })
                funcion_evaluada.append(aproximacion)
            S1 = sum(funcion_evaluada)

            self.pasos.append({ 'titulo':'Calcular S1 con las integrales evaluando los puntos de soporte en posición 0, 3, 4, 7, 10, ...  en los límites de la integral y en la función', 
                            'procedimiento2': procedimiento,
                            'resultado': '\\( S1 = \\sum \\ ' + str(funcion_evaluada) + ' \\ = \\ '+ str(S1) + ' \\)' })

            ##### s2
            funcion_evaluada = []
            procedimiento = []
            puntos_soporteS1 = [soportes[t] for t in range(1,3*particiones,3)]
            for t in puntos_soporteS1:
                aa = float(a.subs(x, t))
                bb = float(b.subs(x, t))
                expresion =  self.exp.subs(x,t)

                integral = integracion_numerica([aa, bb], str(expresion.subs(y, x)) )
                aproximacion = integral.simpson3_8_compuesto(particiones, errores = False)

                procedimiento.append({'titulo':'Calcular evaluando el punto de soporte \\( \\ '+ str(t) +' \\ \\) \\( \\ \\Rightarrow \\int_{'+ str(aa)+ '}^{'+ str(bb)+ '} '+ latex(expresion) +'\\ dy\\ \\) con el método Simpson 3/8 de '+ str(particiones) + ' particiones', 
                                      'procedimiento': integral.pasos,
                                      'resultado': str(aproximacion) 
                                    })
                funcion_evaluada.append(aproximacion)
            S2 = sum(funcion_evaluada)
            

            self.pasos.append({ 'titulo':'Calcular S2 con las integrales evaluando los puntos de soporte en posición 1, 4, 7, 10, ... impares en los límites de la integral y en la función', 
                            'procedimiento2': procedimiento,
                            'resultado': '\\( S2 = \\sum \\ ' + str(funcion_evaluada) + ' \\ = \\ '+ str(S2) + ' \\)' })
            ##### s3
            funcion_evaluada = []
            procedimiento = []
            puntos_soporteS1 = [soportes[t] for t in range(3,3*particiones-1,3)]
            for t in puntos_soporteS1:
                aa = float(a.subs(x, t))
                bb = float(b.subs(x, t))
                expresion =  self.exp.subs(x,t)

                integral = integracion_numerica([aa, bb], str(expresion.subs(y, x)) )
                aproximacion = integral.simpson3_8_compuesto(particiones, errores = False)

                procedimiento.append({'titulo':'Calcular evaluando el punto de soporte \\( \\ '+ str(t) +' \\ \\) \\( \\ \\Rightarrow \\int_{'+ str(aa)+ '}^{'+ str(bb)+ '} '+ latex(expresion) +'\\ dy\\ \\) con el método Simpson 3/8 de '+ str(particiones) + ' particiones', 
                                      'procedimiento': integral.pasos,
                                      'resultado': str(aproximacion) 
                                    })
                funcion_evaluada.append(aproximacion)
            S3 = sum(funcion_evaluada)
            

            self.pasos.append({ 'titulo':'Calcular S3 con las integrales evaluando los puntos de soporte en posición 2, 5, 8, 11, ...  en los límites de la integral y en la función', 
                            'procedimiento2': procedimiento,
                            'resultado': '\\( S3 = \\sum \\ ' + str(funcion_evaluada) + ' \\ = \\ '+ str(S3) + ' \\)' })

            aa1 = float(a.subs(x, c))
            bb1 = float(b.subs(x,c))
            expression1 = self.exp.subs(x,c)

            aa2 = float(a.subs(x, d))
            bb2 = float(b.subs(x, d))
            expression2 = self.exp.subs(x,d)

            procedimiento2 = []
            integral_gc = integracion_numerica([aa1, bb1], str(expression1.subs(y, x)))
            gc = integral_gc.simpson3_8_compuesto(particiones, errores = False)
            procedimiento2.append({'titulo':'Calcular evaluando el punto  \\( \\ '+ str(c) +' \\ \\) \\( \\ \\Rightarrow \\int_{'+ str(aa1)+ '}^{'+ str(bb1)+ '} '+ latex(expression1) +'\\ dy\\ \\) con el método Simpson 3/8 de '+ str(particiones) + ' particiones', 
                                      'procedimiento': integral_gc.pasos,
                                      'resultado': str(gc) 
                                    })
            integral_gd = integracion_numerica([aa2, bb2], str(expression2.subs(y, x)))
            gd = integral_gd.simpson3_8_compuesto(particiones, errores = False)
            procedimiento2.append({'titulo':'Calcular evaluando el punto  \\( \\ '+ str(d) +' \\ \\) \\( \\ \\Rightarrow \\int_{'+ str(aa2)+ '}^{'+ str(bb1)+ '} '+ latex(expression2) +'\\ dy\\ \\) con el método Simpson 3/8 de '+ str(particiones) + ' particiones', 
                                      'procedimiento': integral_gd.pasos,
                                      'resultado': str(gd) 
                                    })

            self.pasos.append({ 'titulo':'Calcular las integrales evaluando los límites y la función con \\( \\ c \\ \\) y \\( \\ d \\ \\) ', 
                            'procedimiento2':procedimiento2 ,
                            'resultado': f'\\(G(c) = {gc} \\ \\ G(d) = {gd} \\)'  })


            self.solucion =  (3*h/8)*(gc + 3*S1 + 3*S2 + 2*S3 + gd ) 
            self.pasos.append({ 'titulo':'Calcular la aproximación con la fórmula', 
                            'procedimiento': '\\( 3 \\cdot h\\cdot\\frac{1}{8} \\cdot (G(c) + G(d) + 3 \\cdot S1 + 3 \\cdot S2 + 2 \\cdot S3) \\)',
                            'resultado': ' \\( \\Rightarrow \\ 3 \\cdot' + str(h) +'\\cdot\\frac{1}{8} \\cdot ('+ str(gc) +' + '+ str(gd) +' + 3 \\cdot '+  str(S1) +'+ 3\\cdot '+ str(S2) + '\\) \n \\( + 2\\cdot '+ str(S3)  +' =  ' + str(self.solucion) + ' \\)'    })
        
            self.metodo = "Simpson 3/8 compuesto doble"
            return N(self.solucion)
            
    ####----- Extrapolación: ------####   
    def romberg(self, n, i=2, metodo = None):
        
        assert n%2 == 0, 'N tiene que ser par'
        
        if i == 2:
            ## O(n^2)
            metodos = {'1':self.trapezoidal_compuesto, 
                       '2':self.simpson1_3_compuesto,
                       '3':self.simpson3_8_compuesto}
            opcion = metodo
            no_trapecios = int(n/2)
                                                                                                          
            self.estimadores = [metodos[opcion](1)] +[metodos[opcion](particiones) for particiones in range(2, no_trapecios*2, 2)]
            
            if n != 2:
                return self.romberg(n, i + 2)
            else:
                return self.estimadores[0]
        
        elif i == n:
            self.solucion = ( 1/(4**((i-2)/2) -1 ) )*((4**((i-2)/2))*self.estimadores[1] - self.estimadores[0])
            self.metodo = "Romberg con " + self.metodo + "y O(h^" + str(n) + ")"
            self.reiniciar_errores()
            return N(self.solucion)
    
        else:
            self.estimadores = [( 1/(4**((i-2)/2) -1 ) )*((4**((i-2)/2))*self.estimadores[pos+1] - self.estimadores[pos]) for pos in range(len(self.estimadores)-1)]
            return self.romberg(n, i + 2)

    ####----- ERRORES: ------####
    def maximo(self, grado, f = None):
        if not f:
            f = self.exp
        x = symbols = 'x'
        g = self.exp
        for _ in range(grado):
            g = diff(g, x)

        try:             
            puntos_criticos = solve(g,x)
        except:
            puntos_criticos = []
        puntos_criticos = [p for p in puntos_criticos if p >= self.a and p<=self.b] + [self.a,self.b]
        
        maximo = float(abs(f.subs(x,puntos_criticos[0])))
        for punto in puntos_criticos:
            valor = float(abs(f.subs(x, N(punto))))
            if valor > maximo:
                maximo = valor
                
        return maximo
    def reiniciar_errores(self):
        self.aproximado = None
        self.estimado = None
        self.cota = None 
        self.total = None
        self.relativo = None
        self.verdadero = None    
    def errores(self, error = None):
        x = symbols('x')
        y = symbols('y')
        if error:
            if error == 'total':
                return self.total
            elif error == 'verdadero':
                return integrate(self.exp, (x, self.a, self.b)) - self.solucion
            elif error == 'relativo ':
                return (1 - self.solucion/integrate(self.exp, (x, self.a, self.b)))*100
            elif error == 'aproximado':
                return self.aproximado            
            elif error == 'estimado':
                return self.estimado
            elif error == 'cota':
                return self.cota
            elif error == 'total':
                return self.total
            else:
                raise ValueError('No existe ese error')
        else: 
            if 'doble' in self.metodo:
                if 'numérico' not in self.metodo:
                    integral_y = integrate(self.exp, (y, self.a, self.b))
                    valor_verdadero = N(integrate (integral_y, (x, self.c, self.d)))
                else:
                    integral_x = integrate(self.exp, (x, self.a, self.b))
                    valor_verdadero = N(integrate (integral_x, (y, self.c, self.d)))
            else:
                valor_verdadero = integrate(self.exp, (x, self.a, self.b))
                
            self.verdadero = valor_verdadero - self.solucion
            self.relativo = (1 - self.solucion/valor_verdadero)*100

            error = {'Total': self.total, 'Verdadero': self.verdadero, 'Relativo': self.relativo, 'Aproximado':self.aproximado, 
                    'Estimado':self.estimado, 'Cota':self.cota}
            return pd.DataFrame(error.values(),index = error.keys(), columns = ['Valor']).reset_index().rename({'index':'Error'}, axis = 1).set_index('Error')