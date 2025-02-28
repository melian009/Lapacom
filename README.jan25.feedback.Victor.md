me he mirado las ecs (1), (2), la fig. 4 y la Tabla.

Para el caso cij=H=0, la solución para el caso simple es N=K(r_i R_i K -
d_i).

Por unidades no puede ser, tiene que ser N= K(r_i R_i -d_i)

d_i está en unidades de y^-1; no tengo ni idea los valores ni unidades
de r_i ni de R_i. R_i es un número entre 0 y 1 y depende de S. Pero en
todo caso en equilibrio es 1. Entonces el único parámetro que nos queda
es r_i.

El resultado es un punto fijo y no una curva. Entonces no entiendo la
curva. Sólo se me ocurre que sea la trayectoria inicializando los
tamaños con S=0. Pero ese caso es muy raro porque los individuos tienen
tamaño cuando nacen pero esta ecuación los hace 0 al principio dela
simulación y no cada vez que un individuo nace.

Además los índices es un lío. En el caso complejo la i se usa pra las
etapas del ciclo pero en el cao simple se una para las especies. Yo
usuaría un índice s para las especies tanto en el caso simple com oel
caso complejo.

Hago una propuesta para el caso simple.

También sugiero quitar la ecuación delos tamaños porque no sirve para
nada. Lo metería como  otro punto la sección 2.3.2. Y diría que la
simulaciones se hacen con R=1.

