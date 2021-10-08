### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 07d14be7-0e98-4745-921c-e0126430d299
begin
	import Pkg
	Pkg.activate(mktempdir())
	
	Pkg.add("PlutoUI")
	using PlutoUI, Plots
	
	Pkg.add("Plots")
	using Plots
	plotly()
	
	Pkg.add("LinearAlgebra") 
	using LinearAlgebra
	
	Pkg.add("NLsolve")
	using NLsolve
	
	Pkg.add("Roots")
	import Roots
end

# ╔═╡ 5d18d960-2461-11ec-2fe8-ed2a60f383c0
md"""
#### Распределение тепла в цилиндре

##### Задание 3, № 10 по книге

Задача симметричная, цилиндр бесконечный, поэтому в каждый момент времени распределение тепла зависит только от радиуса.

Внутри цилиндра расположен концентрический цилиндр, внутри которого происходит выделение тепла с мощностью 
``Q\; \left[ \frac{\text{Вт}}{\text{м}^3} \right]``.
У внутреннего цилиндра другая объёмная теплоёмкость ``C`` и теплопроводность ``k``. 

Уравнение теплопроводности во внутреннем цилиндре:

``\displaystyle
C_i \frac{\partial T}{\partial t} 
= \frac{1}{r} \frac{\partial}{\partial r} \left(
	k_i r \frac{\partial T}{\partial r}
\right) + Q
\qquad \text{(i – значит inner)}``

Во внешнем цилиндре:

``\displaystyle
C_o \frac{\partial T}{\partial t} 
= \frac{1}{r} \frac{\partial}{\partial r} \left(
	k_o r \frac{\partial T}{\partial r}
\right) \phantom{ + q\,}
\qquad \text{(o – outer)}``

Поверхность внешнего цилиндра сообщается со внешней средой (тепловым резервуаром температурой ``T_e``), а также излучает по закону абсолютно чёрного тела. Поэтому граничное условие записывается в виде:

``\displaystyle
\left. 
k_o \frac{\partial T}{\partial r} + 
\alpha \cdot (T - T_e) + \sigma T^4 
\right|_{r = R_2} = 0
``

А ещё в центре цилиндра из соображений симметрии 

``\displaystyle \left.
\frac{\partial T}{\partial r}
\right|_{r=0} = 0 \qquad \forall t``
"""

# ╔═╡ 3ea58b31-7a3a-4cfe-a1c8-3480cc30cc49
md"""
###### Conservation

Согласно уравнению теплопроводности,

``\displaystyle
\frac{\partial T}{\partial t}
\sim
\frac{\partial}{\partial r} \left(
	k r \frac{\partial T}{\partial r}
\right) ``

Следовательно, величина
``\;\displaystyle F(r,t) = k r \frac{\partial T}{\partial r} \;``
должна быть непрерывной по ``r``.

Иначе при дифференциировании возникнет ``\delta``-функция и ``\:\dfrac{\partial T}{\partial t}\:`` (скорость роста температуры) в точке разрыва будет бесконечной.

"""

# ╔═╡ 777760cc-1eac-42eb-ad7b-f79f95f21cc1
md"""
###### Параметры

Пусть ``\alpha = 100 \ \frac{\text{Вт}}{\text{м}^2\text{K}}`` (это типичное значение для воздуха согласно [Wolfram Tutorial on Modeling Heat Transfer](https://reference.wolfram.com/language/PDEModels/tutorial/HeatTransfer/HeatTransfer.html#291584908))

Внутренний цилиндр возьмём из золота, внешний -- из железа.

"""

# ╔═╡ f72cf491-7238-4a61-8a77-34b6503143c1
begin
	
	Tₑ = 0 			# температура окружающей среды, ᵒK
	T₀(r) = 0 		# начальное распределение температур, ᵒK
	
	α = 100. 		# коэф. теплопереноса воздуха, 	Вт / (м²⋅ᵒK)
	σ = 5.67*10^-8 	# постоянная Стефана-Больцмана,	Вт / (м²⋅ᵒK⁴)
	
	kᵢ = 320. 		# коэф. теплопроводности золота, Вт / (мᵒK)
	kₒ = 80. 		# коэф. теплопроводности железа, Вт / (мᵒK)
	
	Cᵢ = 2.49*10^6 	# объёмная теплоёмкость золота, Дж / (м³⋅ᵒK)
	Cₒ = 3.53*10^6 	# объёмная теплоёмкость железа, Дж / (м³⋅ᵒK)
	
	Q = 20000. 		# подводимое тепло, Вт / м³
	
end;

# ╔═╡ a8c158d0-3080-4a1e-8ba3-a9d7d6d4608a
md"""
##### Алгоритм решения

Вдохновляясь работами [(1)](https://www.csc.kth.se/utbildning/kth/kurser/DN2255/ndiff13/Lecture3.pdf) и [(2)](https://krashkov.com/fvm-heat-equation/), применим метод конечных объёмов.

Разобьём цилиндр на цилиндрические слои толщиной ``dr`` и введём интегральное среднее от температуры:

$$u_j(t) = \frac{\displaystyle
\int_{r_{j-1/2}} ^ {r_{j+1/2}}  T(r, t)\, 2\pi r\, dr}{\displaystyle
\int_{r_{j-1/2}} ^ {r_{j+1/2}}  2\pi r\, dr}
= \lambda_j \int\limits_{r_{j-1/2}} ^ {r_{j+1/2}}  T(r, t)\, r\, dr
\;\approx\; T_j,$$

где ``\:\lambda_j = \dfrac{2}{r_{j+1/2}^2 - r_{j-1/2}^2}
\approx \dfrac{1}{r_j h}`` 
только при ``r_j \gg h.``

В методе конечных объёмов вместо значений ``T_j`` нужно вычислять значения ``u_j``.

Рассмотрим область, не включающую в себя границы цилиндров. Коэффициенты ``k`` и ``C`` тогда будут константами.

В силу исходного уравнения, производная по времени от интегрального среднего равна:

``\displaystyle
\frac{d u_j}{dt} 
= \lambda_j \int
	\frac{\partial T}{\partial t} \, r\, dr
= \frac{\lambda_j}{C} \int\limits_{r_{j-1/2}} ^ {r_{j+1/2}}
	\frac{1}{\not{r}} \frac{\partial}{\partial r} \left(
		kr \frac{\partial T}{\partial r}
	\right) \!\! \not{r}\, dr 
= \frac{k\lambda_j}{C} \cdot \left.
	r\, \frac{\partial T}{\partial r}
\right|_{r_{j-1/2}} ^ {r_{j+1/2}}``

Во внутреннем цилиндре есть дополнительное слагаемое:

``\displaystyle
\frac{d u_j}{dt}
= \frac{k\lambda_j}{C} \cdot \left.
	r\, \frac{\partial T}{\partial r}
\right|_{r_{j-1/2}} ^ {r_{j+1/2}} + \frac{Q}{C}
``
"""

# ╔═╡ a4c2549f-2fa9-4186-b2d2-f2d749199789
md"""
###### Аппроксимация пространственных производных

Аппроксимируем поток в полуцелых точках по формуле центральной разности:

``\displaystyle
\left. r \,\frac{\partial T}{\partial r}  \right|_{r_j + 1/2}
\approx r_{j+1/2} \cdot
	\frac{T_{j+1} - T_j}{h}
\approx r_{j+1/2} \cdot
	\frac{u_{j+1} - u_j}{h}
``

Подставляя в полученную выше формулу, получим:

``\displaystyle
\frac{d u_j}{dt}
= \frac{k\lambda_j}{C} \cdot \left.
	r\, \frac{\partial T}{\partial r}
\right|_{r_{j-1/2}} ^ {r_{j+1/2}} 
\approx \dfrac{k\lambda_j}{Ch} \Bigl(
	r_{j-1/2} \cdot u_{j-1}
    - (r_{j-1/2} + r_{j+1/2}) \cdot u_j
	+ r_{j+1/2} \cdot u_{j+1}
\Bigr)``

``
\phantom{\dfrac{d u_j}{dt}}
= \dfrac{k\lambda_j}{Ch} \Bigl(
	r_{j-1/2} \cdot u_{j-1}
    - 2 r_j \cdot u_j
	+ r_{j+1/2} \cdot u_{j+1}
\Bigr)``

Согласно работе [(1)](https://www.csc.kth.se/utbildning/kth/kurser/DN2255/ndiff13/Lecture3.pdf) порядок аппроксимации у аналогичной схемы ``\;-\; \displaystyle O(h^2)``.

"""

# ╔═╡ fe6815e3-9300-441c-b628-db76c0ac1dc4
md"""
###### Учёт непрерывности потока

Расположим один из узлов на границе между цилиндрами. Обозначим его номер как ``\xi``.

Так как величина ``\:kr\,\dfrac{\partial T}{\partial r}\:`` должна быть непрерывной,

$$\left. k_i \dfrac{\partial T}{\partial r} \right|_{r = r_\xi-0}
=
  \left. k_o \dfrac{\partial T}{\partial r} \right|_{r = r_\xi+0}$$

Используя формулы дифференциирования вперёд/назад:

$$k_i \cdot \bigl( u_{\xi-2} - 4u_{\xi-1} + 3u_\xi \bigr)
= k_o \cdot \bigl(-3u_\xi + 4u_{\xi+1} - u_{\xi+2} \bigr)$$

$$u_\xi = \frac{1}{3(k_o + k_i)} \bigl(
	- k_i \cdot u_{\xi-2} + 4k_i \cdot u_{\xi-1}
	+4k_o \cdot u_{\xi+1} -  k_o \cdot u-{\xi+2}
\bigr) \tag{1}$$

"""

# ╔═╡ 08c155dd-f25d-4211-b545-bb028aedc42a
md"###### Аппроксимация граничных условий"

# ╔═╡ 1a4cf5e3-ca38-4777-90b0-30db51d7a9fc
md"""
1). ``\quad \displaystyle \left.
\frac{\partial u}{\partial r}
\right|_{r=0} = 0 ``

Разместим первый узел в точке ``\displaystyle r_1 = h/2``, а также добавим фиктивный нулевой узел в точке ``\displaystyle r_0 = -h/2``.

Разностная схема в первой точке запишется в виде:

``\displaystyle \frac{d u_1}{dt}
\approx \frac{k_i\lambda_1}{C_i h} \Bigl(
	\underbrace{r_{1-1/2}}_{0} \cdot u_0
	- \underbrace{2r_1}_{h} \cdot u_1 
	+ \underbrace{r_{1+1/2}}_{h} \cdot u_2
\Bigr) + \frac{Q}{C_i}
\qquad\quad
\biggl[ \lambda_1 = \frac{2}{h^2} \biggr]
``

``\displaystyle \phantom{\frac{d u_1}{dt}}
= \frac{2k_i}{C_i h^2} \Bigl(
	-u_1 + u_2
\Bigr) +\frac{Q}{C_i}
``

Напомню, что индекс ``i`` означает внутренний цилиндр.

Значение в фиктивной нулевой точке оказалось неважным. Левое граничное условие аппроксимировать не удалось.

Надеемся, что оно как-нибудь выполнится само собой.
Если всё же не выполнится, то можно будет использовать формулу дифференциирования вперёд.
"""

# ╔═╡ 5246439f-68b6-4253-9288-561371ca55ba
md"""
2). ``\quad\displaystyle
\left. 
k_o \frac{\partial T}{\partial r} + 
\alpha \cdot (T - T_e) + \sigma T^4 
\right|_{r = R_2} = 0
``

Расположим узел последний узел сетки ровно на поверхности внешнего цилиндра.

Введём фиктивный узел ``r_{N+1} = R_2 + h``. Тогда по формуле центральной разности:

``k_o \dfrac{u_{N+1} - u_{N-1}}{2h} + 
\alpha \cdot (u_N - T_e) + \sigma u_N^4 = 0``

``
u_{N+1} = u_{N-1} - \dfrac{2h}{k_o} \Bigl(\alpha \cdot (u_N - T_e) + \sigma u_N^4 \Bigr)
``

В случае явной схемы можно вычислить ``\:u_{N+1},\:`` а затем вычислить ``\:\hat u_N\:`` по общему шаблону.
"""

# ╔═╡ ab252169-1ce2-49d7-9dc1-838b656f3128
md"""
###### Аппроксимация временно́й производной

| Метод          | Устойчивость                                       | Монотонность |
|----------------|:--------------------------------------------------:|:------------:|
|Явный Эйлера    |``\frac{k}{C} \cdot \frac{\tau}{h^2} \leqslant \frac{1}{2}``|``\frac{k}{C} \cdot \frac{\tau}{h^2} < \frac{1}{2}``|
|Неявный Эйлера  |      ``\forall``       |      ``\forall``       |
|Кранка-Николсона|      ``\forall``       |``\frac{k}{C} \cdot \frac{\tau}{h^2} \lesssim \frac{1}{2}``|

Выберу неявный метод Эйлера, поскольку он безусловно устойчивый и монотонный.


Шаблон:

```
* - * - *
    |
    *
```

"""

# ╔═╡ 3131b45a-21e1-4bb7-b0e6-136675746b6c
md"""

###### Составим систему уравнений

Точкам ``\:1\ldots\xi-1,\quad \xi+1\ldots N\:`` (т.е. кроме точки разрыва) соответствуют уравнения

$$\begin{gathered}\dfrac{\hat{u} - u}{\tau}
= \frac{k\lambda_j}{Ch} \Bigl(
	r_{j-1/2} \cdot \hat{u}_{j-1}
    - 2 r_j \cdot \hat{u}_j
	+ r_{j+1/2} \cdot \hat{u}_{j+1}
\Bigr) + \frac{Q}{C}\\
\Leftrightarrow\\[0.2em]
-\frac{k\lambda_j\tau}{Ch} r_{j-1/2} \cdot \hat{u}_{j-1}
+ \Bigl(1 + 2 \frac{k\lambda_j\tau}{Ch} r_j \Bigr) \cdot \hat{u}_j
-\frac{k\lambda_j\tau}{Ch} r_{j+1/2} \cdot \hat{u}_{j+1}
= u_j + \frac{Q}{C}\tau
\end{gathered}$$

Точке разрыва соответствует уравнение ``(1)``:

$$\frac{k_i}{3(k_o + k_i)} \hat{u}_{\xi-2} 
- \frac{4}{3}\frac{k_i}{k_o + k_i} \hat{u}_{\xi-1}
+ \hat{u}_\xi
- \frac{4}{3}\frac{k_o}{k_o + k_i} \hat{u}_{\xi+1} 
+ \frac{k_o}{3(k_o + k_i)} \hat{u}_{\xi+2} 
= 0$$

Последней, фиктивной точке соответствует уравнение

$$-u_{N-1} + \frac{2h}{k_o}\bigl(
	\alpha u_N + \sigma u_N^4
\bigr) 
+u_{N+1} = \frac{2h}{k_o}\alpha T_e$$

Вместе эти уравнения образуют нелинейную систему, которую мы будем решать с помощью пакета NLsolve.jl методом Ньютона.
"""

# ╔═╡ bf661b18-5915-4658-89cd-2764a3117c54
md"""Заметим, что матрица Якоби этой системы практически трёхдиагональная, мешается только последняя строка и строка номер ``\xi``.

Исправим это с помощью линейных преобразований.

Рассмотрим уравнения под номером ``\;\xi-1,\quad \xi,\quad \xi+1``: 

\
$$\left[\begin{array}{ccccc|c}
-\frac{k_i \lambda_{\xi-1}\tau}{C_i h}r_{\xi-3/2} &
1 + 2\frac{k_i \lambda_{\xi-1}\tau}{C_i h}r_{\xi-1} &
-\frac{k_i \lambda_{\xi-1}\tau}{C_i h}r_{\xi-1/2} &
⋅ & ⋅ &
u_{\xi-1} + \frac{Q}{C}\tau
\\
\frac{k_i}{3(k_o + k_i)} &
-\frac{4}{3}\frac{k_i}{k_o+k_i} &
1 &
-\frac{4}{3}\frac{k_o}{k_o+k_i} &
\frac{k_o}{3(k_o + k_i)} &
0
\\
⋅ & ⋅ &
-\frac{k_o \lambda_{\xi+1}\tau}{C_o h}r_{\xi+1/2} &
1 + 2\frac{k_o \lambda_{\xi+1}\tau}{C_o h}r_{\xi+1} &
-\frac{k_o \lambda_{\xi+1}\tau}{C_o h}r_{\xi+3/2} &
u_{\xi+1}
\end{array}\right]$$

$$\sim$$

Прибавим к средней строчке первую и третью строчки так, чтобы часть коэффициентов занулилась и матрица стала трёхдиагональной:

$$\left[\begin{array}{ccccc|c}
-\frac{k_i \lambda_{\xi-1}\tau}{C_i h}r_{\xi-3/2} &
1 + 2\frac{k_i \lambda_{\xi-1}\tau}{C_i h}r_{\xi-1} &
-\frac{k_i \lambda_{\xi-1}\tau}{C_i h}r_{\xi-1/2} &
⋅ & ⋅ &
u_{\xi-1} + \frac{Q}{C}\tau
\\
⋅ &
\frac{2}{3}\frac{k_i}{k_o+k_i}\bigl(\frac{r_{\xi-1}}{r_{\xi-3/2}} - 2\bigr) + \frac{C_ih}{3(k_o+k_i)\lambda_{\xi-1}\tau r_{\xi-3/2}} &
1 - \frac{1}{3(k_o+k_i)}\bigl(k_i\frac{r_{\xi-1/2}}{r_{\xi-3/2}} + k_o\frac{r_{\xi+1/2}}{r_{\xi+3/2}}\bigr) &
\frac{2}{3}\frac{k_o}{k_o+k_i}\bigl(\frac{r_{\xi+1}}{r_{\xi+3/2}}-2\bigr) + \frac{C_oh}{3(k_o+k_i)\lambda_{\xi+1}\tau r_{\xi+3/2}}  &
⋅ &
\frac{h}{3\tau(k_o+k_i)}\left(\frac{C_i}{\lambda_{\xi-1} r_{\xi-3/2}} \bigl(u_{\xi-1} + \frac{Q}{C_i}\tau \bigr) +
\frac{C_o}{\lambda_{\xi+1} r_{\xi+3/2}} u_{\xi+1} \right)
\\
⋅ & ⋅ &
-\frac{k_o \lambda_{\xi+1}\tau}{C_o h}r_{\xi+1/2} &
1 + 2\frac{k_o \lambda_{\xi+1}\tau}{C_o h}r_{\xi+1} &
-\frac{k_o \lambda_{\xi+1}\tau}{C_o h}r_{\xi+3/2} &
u_{\xi+1}
\end{array}\right]$$

"""

# ╔═╡ 416513b5-ad52-46f7-849e-52fa9cc9f628
md"""
Рассмотрим уравнения под номером ``\:\;N,\quad N+1``:

$$\left\{\begin{array}{cccl}
-\frac{k_o\lambda_N\tau}{C_oh} r_{N-1/2} \cdot \hat{u}_{N-1}
& + \Bigl(1 + 2 \frac{k\lambda_N\tau}{C_oh} r_N \Bigr) \hat{u}_N
& -\frac{k\lambda_N\tau}{C_oh} r_{N+1/2} \cdot \hat{u}_{N+1}
& = u_N
\\
-\hat{u}_{N-1} 
& + \frac{2h}{k_o}\bigl(\alpha u_N + \sigma \hat{u}_N^4\bigr) 
& +\hat{u}_{N+1} 
& = \frac{2h}{k_o}\alpha T_e
\end{array}\right.$$

$$\Leftrightarrow$$

$$\left\{\begin{array}{cccl}
-\frac{k_o\lambda_N\tau}{C_oh} r_{N-1/2} \cdot \hat{u}_{N-1}
& + \Bigl(1 + 2 \frac{k\lambda_N\tau}{C_oh} r_N \Bigr) \hat{u}_N
& -\frac{k\lambda_N\tau}{C_oh} r_{N+1/2} \cdot \hat{u}_{N+1}
& = u_N
\\

& \frac{2h}{k_o}\bigl(\alpha u_N + \sigma \hat{u}_N^4\bigr)
-\frac{C_oh}{k_o\lambda_N\tau r_{N-1/2}} \hat{u}_N
- 2 \frac{r_N}{r_{N-1/2}} u_N
& +u_{N+1} + \frac{r_{N+1/2}}{r_{N-1/2}} \hat{u}_{N+1} 
& = -\frac{C_oh}{k_o\lambda_N\tau r_{N-1/2}} \hat{u}_N 
+\frac{2h}{k_o}\alpha T_e 
\end{array}\right.$$

Якобиан системы совпадает с матрицей коэффициентов во всех строчках, кроме последней, потому что соответствующие уравнения линейные. Последняя же строчка якобиана равна:

$$\begin{bmatrix}
0, & \ldots & 0, & 
\frac{2h\alpha}{k_o} + 3\sigma u_N^4 - \frac{C_oh}{k_o\lambda_N\tau r_{N-1/2}} - 2 \frac{r_N}{r_{N-1/2}}, &
1 + \frac{r_{N+1/2}}{r_{N-1/2}} 
\end{bmatrix}$$
"""

# ╔═╡ 0ecaeec8-b990-4760-add6-9e3701a98b35
begin
	Nx= 1000 		# Число узлов, не считая фиктивного узла № (Nx+1)
	ξ = 500 		# Номер узла, через который проходит граница между цилиндрами
	h = 0.0001 		# шаг по радиусу в метрах
	R₁= -h/2 + ξ*h 	# Радиус внутреннего цилиндра (первый узел находится в r = h/2)
	R₂= -h/2 + Nx*h # Радиус внешнего цилиндра
	
	τ = 100. 			# шаг по времени в секундах
	
	T = 20000 			# Время моделирования в секундах
	Nt=Int(ceil(T/τ)) 	# Число точек по времени
	T = τ * Nt 			# Уточнённое время моделирования (если не делится нацело)
end;

# ╔═╡ f2b172f6-10af-4529-a00e-742847909a32
(τ, Nt, h)

# ╔═╡ 21a681fe-0998-4f77-b497-f38c7acf5bb5
Nt * (Nx+1) * 8 / 1024 / 1024

# ╔═╡ ebed84dd-f040-47a9-ad9d-8ee200f1ec3f
md"#### Решение. Код"

# ╔═╡ 09a24c87-db8a-453f-babf-9f67ccd265f9
# Функция, возвращающая координаты узлов сетки (неважно, целых или полуцелых)

r(j) = h * (j - 1/2);

# ╔═╡ 0a1a2660-5b98-44dd-952f-bfaa60032c62
params = (Nt, Nx, ξ, τ, h, kᵢ, kₒ, Cᵢ, Cₒ, Q, T₀, r);

# ╔═╡ 40428de5-ecc7-4152-b058-0cd23d9f66b2
function calcU(params)
	
	Nt, Nx, ξ, τ, h, kᵢ, kₒ, Cᵢ, Cₒ, Q, T₀, r = params
	
# Инициализация массивов
	
	# Здесь будет храниться решение
	u = Array{Float64}(undef, Nx+1, Nt+1)
	
	# ~Обратные площади колец, более точно см. начало описания метода конечных объёмов
	λ = [2 / (r(j+1/2)^2 - r(j-1/2)^2) for j = 1 : Nx]
	
	# Первый временной слой:
	u[:, 1] = (T₀ ∘ r).(1 : Nx+1)
	
	
	# Вектор правых частей 
	# (будем с ним решать систему уравнений на каждом временно́м шаге)
	rhs = Array{Float64}(undef, Nx+1)
	

	
# Создадим theMatrix - постоянную часть Якобиана:
	#
	# Диагональ под главной диагональю
	dl = Array{Float64}(undef, Nx)
	dl[1:ξ-2]  = [-kᵢ*λ[j]*τ/Cᵢ/h*r(j-1/2) for j ∈ 2 : ξ-1]
	dl[ξ-1] = 2/3*kᵢ/(kᵢ+kₒ)*(r(ξ-1)/r(ξ-3/2)-2) + Cᵢ*h / 3(kₒ+kᵢ)/λ[ξ-1]/τ/r(ξ-3/2)
	dl[ξ:Nx-1] = [-kₒ*λ[j]*τ/Cₒ/h*r(j-1/2) for j ∈ ξ+1 : Nx]
	#
	# Главная диагональ
	d = Array{Float64}(undef, Nx+1)
	d[1 : ξ-1]  = [1 + 2kᵢ*λ[j]/Cᵢ * τ/h * r(j) for j ∈ 1 : ξ-1]
	d[ξ] = 1 - 1/3(kᵢ+kₒ)*(kᵢ*r(ξ-1/2)/r(ξ-3/2) + kₒ*r(ξ+1/2)/r(ξ+3/2))
	d[ξ+1 : Nx] = [1 + 2kₒ*λ[j]/Cₒ * τ/h * r(j) for j ∈ ξ+1 : Nx]
	d[Nx+1] = 1 + r(Nx+1/2)/r(Nx-1/2)
	#
	# Диагональ над главной
	du = Array{Float64}(undef, Nx)
	du[1:ξ-1]  = [-kᵢ*λ[j]*τ/Cᵢ/h*r(j+1/2) for j ∈ 1 : ξ-1]
	du[ξ] = 2/3*kₒ/(kᵢ+kₒ)*(r(ξ+1)/r(ξ+3/2)-2) + Cₒ*h / 3(kₒ+kᵢ)/λ[ξ+1]/τ/r(ξ+3/2)
	du[ξ+1:Nx] = [-kₒ*λ[j]*τ/Cₒ/h*r(j+1/2) for j ∈ ξ+1 : Nx]
	#
	theMatrix = Tridiagonal(dl, d, du)
	
	
	
# Модулю NLsolve.jl нужны функции f! и j!, которые вычисляют невязку и Якобиан
	
	# Невязка Aû-rhs:
	# Первые Nx элементов можно просто получить умножением на матрицу.
	# А последний элемент задаётся нелинейным уравнением (краевым условием, где σT⁴),
	# поэтому этот элемент перезаписываем
	#
	function f!(F, û)
		F[:] = theMatrix * û
		F[Nx+1] = 2h/kₒ*(α*û[Nx] + σ*û[Nx]^4) - 
			Cₒ*h / kₒ/λ[Nx]/τ/r(Nx-1/2) * û[Nx] - 2r(Nx)/r(Nx-1/2)*û[Nx] +
			+û[Nx+1] + r(Nx+1/2)/r(Nx-1/2)*û[Nx+1]
		F[:] -= rhs # rhs уже будет вычислена 
					# она не зависит от итераций метода Ньютона (û), 
					# зависит только от предыдущего временного слоя (ù)
	end

	
	# Якобиан:
	# Представляет собой практически постоянную матрицу
	# Нужно обновить единственный элемент
	function j!(J, û)
		J[Nx+1, Nx] = h*α/kₒ + 3σ*û[Nx]^4 - 
					  Cₒ*h / kₒ/λ[Nx]/τ/r(Nx-1/2) - 
					  2r(Nx)/r(Nx-1/2)
	end
	
	
	# See https://github.com/JuliaNLSolvers/NLsolve.jl#if-the-jacobian-is-sparse
	initial_x = Array{Float64}(undef, Nx+1)
	initial_f = Array{Float64}(undef, Nx+1)
	df = OnceDifferentiable(f!, j!, initial_x, initial_f, theMatrix)

	
# The Cycle
	
	for t ∈ 1 : Nt
		
		ù = @view u[:, t] # Старый временной слой
		
		# Вычислим вектор правых частей
		rhs[:] = ù[:]
		rhs[1 : ξ-1] .+= τ*Q/Cᵢ
		rhs[ξ] = h/τ/3(kₒ+kᵢ) * (Cᵢ/λ[ξ-1]/r(ξ-3/2) * (ù[ξ-1] + τ*Q/Cᵢ)
		   + Cₒ/λ[ξ+1]/r(ξ+3/2) * ù[ξ+1])
		rhs[Nx+1] = - Cₒ*h / kₒ/λ[Nx]/τ/r(Nx-1/2) * ù[Nx] +
				+ 2h/kₒ*α*Tₑ
		
		# Решим систему уравнений -→ следующий временно́й слой:
		u[:, t+1] .= nlsolve(df, ù, method = :newton).zero

	end
	
	return u
end;

# ╔═╡ 45bc164a-5693-40c0-bc49-cbf952159b02
u = calcU(params);

# ╔═╡ bf4c3994-cabf-47cb-869a-bf9fd4486b72
@bind t Slider(1:Nt+1)

# ╔═╡ a6cd33ef-3880-427b-b3d9-4e3e008cf7cd
md"""
#### Аналитическое решение при ``\;t \to \infty``

В установившемся режиме ``\dfrac{\partial T}{\partial t} = 0``.

Уравнение становится следующим:

``
\dfrac{1}{r} \dfrac{\partial}{\partial r} \left(
	k r \dfrac{\partial T}{\partial r}
\right) = - Q
``

`` k r \dfrac{\partial T}{\partial r}
= -\dfrac{Q}{2} r^2 + A ``

``
T(r) = -\dfrac{1}{4} \dfrac{Q}{k}\, r^2 + \dfrac{A}{k}\, \ln(r) + B ``,

где ``A`` и ``B`` -- произвольные константы.

Во внутреннем цилиндре решение лучше записать в виде

``T_{in}(r) = -\dfrac{1}{4} \dfrac{Q}{k_i}\, r^2 + B_i``

А во внешнем

``T_{out}(r) = \dfrac{A}{k_o}\, \ln\Bigl(\dfrac{r}{R_2} \Bigr) + B_o``

Приравнивая потоки на границе раздела, получим ``\:A = -\dfrac{1}{2} QR^2_1``

``B_o`` найдём численно из правого граничного условия:

``\dfrac{A}{R_2} + \alpha \cdot (B_o - T_e) + \sigma B_o^4 = 0 ``

Приравняв значения температуры на границе, получим:

``\:B_i = B_o + QR_1^2 \cdot \left( \dfrac{1}{4k_i} - \dfrac{1}{2k_o}\, \ln\Bigl(\dfrac{R_1}{R_2} \Bigr) \right) ``

"""

# ╔═╡ 03234a65-63d6-4a93-aa7b-e28f9ff1c996
begin
	A = -1/2 * Q * R₁^2
	Bₒ = Roots.find_zero(T -> A/R₂ + α*(T-Tₑ) + σ*T^4, Tₑ)
	Bᵢ = Bₒ + Q * R₁^2 * (1/4kᵢ - 1/2kₒ * log(R₁/R₂))
end;

# ╔═╡ 35a6cd9c-2a96-4f9d-a7a0-00c06152c87c
anal(r) =
	if r < R₁
		-1/4 * Q/kᵢ * r^2 + Bᵢ
	else
		A/kₒ * log(r/R₂) + Bₒ
	end

# ╔═╡ c02c0478-99c6-44a4-ae02-c02b3830973a
xlist = r.(1:Nx);

# ╔═╡ 71c26f9d-b2f5-491e-a92f-7d5cb0d5db36
begin
	plot(xlist, anal.(xlist),
		title = "t = $(Int(round(((t-1)*τ)))) с",
		label = "Аналитическое при t = ∞",
		xlabel = "r, м"
	)
	plot!(xlist, u[1 : end-1, t],
		label = "Численное при заданном t",
		#ylims = (270.1, 270.63),
	)
end

# ╔═╡ Cell order:
# ╠═07d14be7-0e98-4745-921c-e0126430d299
# ╟─5d18d960-2461-11ec-2fe8-ed2a60f383c0
# ╟─3ea58b31-7a3a-4cfe-a1c8-3480cc30cc49
# ╟─777760cc-1eac-42eb-ad7b-f79f95f21cc1
# ╠═f72cf491-7238-4a61-8a77-34b6503143c1
# ╟─a8c158d0-3080-4a1e-8ba3-a9d7d6d4608a
# ╟─a4c2549f-2fa9-4186-b2d2-f2d749199789
# ╟─fe6815e3-9300-441c-b628-db76c0ac1dc4
# ╟─08c155dd-f25d-4211-b545-bb028aedc42a
# ╟─1a4cf5e3-ca38-4777-90b0-30db51d7a9fc
# ╟─5246439f-68b6-4253-9288-561371ca55ba
# ╟─ab252169-1ce2-49d7-9dc1-838b656f3128
# ╟─3131b45a-21e1-4bb7-b0e6-136675746b6c
# ╟─bf661b18-5915-4658-89cd-2764a3117c54
# ╟─416513b5-ad52-46f7-849e-52fa9cc9f628
# ╠═0ecaeec8-b990-4760-add6-9e3701a98b35
# ╠═f2b172f6-10af-4529-a00e-742847909a32
# ╠═21a681fe-0998-4f77-b497-f38c7acf5bb5
# ╟─ebed84dd-f040-47a9-ad9d-8ee200f1ec3f
# ╠═09a24c87-db8a-453f-babf-9f67ccd265f9
# ╠═0a1a2660-5b98-44dd-952f-bfaa60032c62
# ╠═40428de5-ecc7-4152-b058-0cd23d9f66b2
# ╠═45bc164a-5693-40c0-bc49-cbf952159b02
# ╟─bf4c3994-cabf-47cb-869a-bf9fd4486b72
# ╠═71c26f9d-b2f5-491e-a92f-7d5cb0d5db36
# ╟─a6cd33ef-3880-427b-b3d9-4e3e008cf7cd
# ╠═03234a65-63d6-4a93-aa7b-e28f9ff1c996
# ╠═35a6cd9c-2a96-4f9d-a7a0-00c06152c87c
# ╠═c02c0478-99c6-44a4-ae02-c02b3830973a
