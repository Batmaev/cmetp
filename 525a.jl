### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ bff05cbb-0331-47c5-be77-d8c3427eb8d6
begin
	import Pkg
	Pkg.activate(mktempdir())

	Pkg.add("LinearAlgebra")
	using LinearAlgebra
	Pkg.add("Plots")
	using Plots
	plotly()
	Pkg.add("OffsetArrays")
	using OffsetArrays
end;

# ╔═╡ 5f3959b0-865c-11ec-13e2-e72527971dcb
md"""
## Задача 5, № 25 по книге

Стационарное уравнение Шрёдингера в одномерной прямоугольной потенциальной яме:
```math
- \frac12 u''(x) + (V(x) - E)\, u(x) = 0 \tag{1}
```
```math
V(x) = 
\begin{cases}
-V_0, & |x| < 1 \\
0, & |x| > 1
\end{cases}
```
Требуется численно найти собственные функции ``u(x)`` и собственные числа ``\gamma``.

Ограничимся связанными состояниями ``(\gamma < 0)``. Тогда вне ямы волновая функция будет экспоненциально убывать, можно будет ограничиться конечным промежутком
```math
x \in [-a,\: a], \;\; a \gg 1 
```
и ввести граничные условия
```math
u(a) = 0 = u(-a)
```
"""

# ╔═╡ b7667130-f054-4962-8937-8bd50b0aab4c
md"""
##### Метод конечных элементов """

# ╔═╡ d2e5fed7-1f31-4dac-8a25-1166ad1aa2d0
md"""
Умножим уравнение (1) на некоторую функцию ``f_k(x)``, зануляющуюся на концах отрезка ``[-a,\: a]``, и проинтегрируем:
```math
\begin{align*}
-\frac12 \int_{-a}^a u'' f_k\, dx + \int_{-a}^a (V(x) - E)uf_k\, dx & = 0 \\
-\cancel{u' f_k \biggr|_{-a}^{a}} + \int_{-a}^a \Bigl(u' f'_k  + 2(V-E)uf_k\Bigr)dx &= 0 \tag{2}
\end{align*}
```
"""

# ╔═╡ 00b7ca39-7a84-4a44-ac09-02e587b1aae5
md"""
Введём сетку 
```math 
x_j = h \cdot j, \quad j = -N..N, \quad h = \frac{a}{N+1} \quad \text{(нет точек ±a)}
```
и базисные функции 
```math
f_j(x) = 
\begin{cases}
1 - \frac{|x - x_j|}{h}, & |x - x_j| < h \\
0, & |x - x_j| > h
\end{cases}
```

Будем искать решение в виде 
```math
u(x) = \sum_j u_j f_j(x) 
\tag{3}
```
Тогда автоматически будет выполнено, что
```math
u(x_j) = u_j, \qquad u(a) = 0 = u(-a)
```
Кроме того будет выполнено
```math
u''(x) \equiv 0
```
что противоречит исходному уравнению (1), но мы закроем на это глаза.
"""

# ╔═╡ 6b1e8f85-e38f-435f-b865-0e505aa16bd2
md"""
Подставим анзац (3) в проинтегрированное уравнение (2):
```math
\sum_j u_j \int_{-a}^a \biggl( f'_j f'_k + 2(V - E) f_j f_k \biggr)dx = 0
\tag{4}
```
Обозначим
```math
A_{jk} = \int_{-a}^a \bigl( f'_j f'_k + 2V(x)f_j f_k \bigr)dx \qquad 
B_{jk} = 2\int_{-a}^a f_j f_k \, dx
```
Тогда последнее уравнение перепишется в виде
```math
\hat A \vec u = E \cdot \hat B \vec u \tag{5}
```
Это *обобщённая задача поиска собственных значений,* решение которой уже реализовано в стандартной библиотеке языка Julia.
"""

# ╔═╡ 166fd752-34c1-41b9-8032-26c86a7ba725
md"""
##### Осталось найти матрицы ``A`` и ``B``.
"""

# ╔═╡ af310304-aa61-49ef-a1db-90c8abd89dbf
md"""
Диагональные элементы B:

``\displaystyle
B_{jj} = 2\int_{-h}^h \left(1 - \frac{|t|}{h} \right)^2 dt = \frac43 h
``

Элементы на двух соседних диагоналях:

``\displaystyle
B_{j,j-1} = B_{j,j+1} = 2\int_0^h \Bigl(1 - \frac{t}{h} \Bigr) \cdot \frac{t}{h}\, dt = \frac13 h
``

Остальные элементы равны ``0``.
"""

# ╔═╡ 802591f2-c124-4663-9c30-b0ebc0c0c2ce
md"""
###### Матрица ``A``
"""

# ╔═╡ 85e2e1dd-c618-43d8-bedd-022ebc86ee24
md"""

Пусть ``\, h = 1/M,\ M \in \mathbb{N},\,`` тогда границы ямы расположены в узлах сетки под номером ``\pm M``.

Диагональные элементы матрицы ``A``:

``\displaystyle
A_{jj} = 
\begin{cases}
\displaystyle\int_{-h}^h \left[\left(\pm \frac1h\right)^2 - 2V_0 \left(1-\frac{|t|}{h} \right)^2 \right]
dt = \frac2h - \frac43 h V_0 & \qquad |j| < M \\[0.5em] 
\displaystyle\int_{-h}^h \left(\pm \frac1h\right)^2 dt = \frac2h & \qquad |j| > M \\[0.5em] 
\displaystyle\int_{-h}^0 (\dots)\, dt + \displaystyle\int_{0}^h (\dots)\, dt = \frac2h - \frac23 h V_0 & \qquad |j| = M
\end{cases}
``

"""

# ╔═╡ a9152c5a-30c1-4120-a8ff-523b10a8b645
md"""
Диагональ над главной:

``\displaystyle
A_{j,j+1} = 
\begin{cases}
\displaystyle\int_0^h \left[- \frac{1}{h^2} - 2V_0 \Bigl(1 - \frac{t}{h} \Bigr) \cdot \frac{t}{h} \right] dt
= - \frac1h - \frac13 h V_0  & \qquad 
-M < j+1 \leqslant M  \\[0.5em] 
\displaystyle\int_{0}^h -\frac1{h^2} \, dt = - \frac1h & \qquad j+1 \leqslant -M \;\;||\;\; j \geqslant M 
\end{cases}
``

Диагональ под главной отдельно вычислять не будем, т.к. матрицы ``A`` и ``B`` симметричные.

"""

# ╔═╡ 30efff26-9f0e-4ef3-8c3d-668a44e8c849
V₀ = 16
# Если задать слишком маленькое значение, 
# может понадобиться расширить область интегрирования
# (Вне ямы решение убывает примерно как e^-x√2V₀.
# Требуется, чтобы оно занулялось на концах)

# ╔═╡ 9806ac50-5f3e-470e-a5ef-6b33b285712d
# Сетка:

begin
	M = 32						# Целое число такое, что h = 1/M
	h = 1/M

	a_min = 8 					# Число такое, что a ⩾ a_min
	a = h * ceil(a_min / h) 	# При |x| > a решение будет считаться равным 0

	N = round(Int, a/h) - 1 	# Точки нумеруются от -N до N
end

# ╔═╡ 3d9ed0c9-5f9e-4514-9593-8637c466177c
params = (; V₀, a_min, M)

# ╔═╡ e761b552-bdce-4abe-a6aa-f12a1302fc3e
function shrodingerSolve(; V₀, a_min, M)
	# Сетка:
	h = 1/M
	a = h * ceil(a_min / h)

	N = round(Int, a/h) - 1
	
	# Матрица B
	b_main_di   = [4h/3 for _ in 1 : 2N+1]
	b_second_di = [ h/3 for _ in 1 : 2N  ]
	B = SymTridiagonal(b_main_di, b_second_di)

	# Матрица A
	a_main_di = OffsetArray( Vector{Float64}(undef, 2N+1),  -N:N )
	a_main_di[-M+1 : M-1] .= 2/h - 4/3 * h * V₀
	a_main_di[-M]          = 2/h - 2/3 * h * V₀
	a_main_di[+M] 		   = 2/h - 2/3 * h * V₀
	a_main_di[-N : -M-1]  .= 2/h
	a_main_di[M+1 : N]    .= 2/h

	a_second_di = OffsetArray( Vector{Float64}(undef, 2N),  -N+1 : N )
	a_second_di[begin : -M ] .= -1/h
	a_second_di[-M+1  :  M ] .= -1/h - 1/3 * h * V₀
	a_second_di[ M+1  : end] .= -1/h

	A = SymTridiagonal(
		OffsetArrays.no_offset_view(a_main_di),
		OffsetArrays.no_offset_view(a_second_di)
	)

	# Решение
	return eigen!(Array(A), Array(B))
end

# ╔═╡ 23a1d88c-756a-46e3-9ea4-d70a9260ebc7
E, u = shrodingerSolve(; params...)

# ╔═╡ 8b433c95-8b23-4c05-881a-fc53249472b6
md"""
Получили $(2N + 1) решений. Мы строили решение в предположении, что энергия отрицательна и волновые функции зануляются на бесконечности. Решения с положительной энергией нужно поэтому отбросить.
"""

# ╔═╡ bcbc3be4-6414-4d9d-ad12-e28edc25844a
begin
	n = 4 # Номер уровня
	
	plot(h.*(-N : N),
		u[:, n],
		xlims = (-3.5, 3.5),
		title = "E<sub>$n</sub> = $(E[n])",
		ylabel = "ψ (у.е.)",
		xlabel = "Координата (у.е.)<br>Яма простирается от -1 до 1"
	)
end

# ╔═╡ d03d4e6e-4f14-41e4-a37c-a345783ea479
md"""
##### Сравнение с теоретическим решением

Задачу практически полностью можно решить аналитически (см. Бегларян Лилит, ТПФ, 6 семестр, квантовая механика, зад. 2а). 

Существуют 2 семейства решений: чётное и нечётное.

``\displaystyle
u_\text{even}(x) = 
\begin{cases}
C_1 \cos(kx), 			& |x| < a \\ 
C_2 e^{-\varkappa |x|}, & |x| \geqslant a
\end{cases}
\qquad\quad
k : \ \sqrt{\frac{2mV_0 a^2}{\hbar^2} - (ka)^2} = (ka) \tan(ka)
``

``\displaystyle \: 
u_\text{odd}(x) = 
\begin{cases}
C_1 \sin(kx), 	       & |x| < a \\ 
C_2 e^{\varkappa x},   & x < -a \\ 
-C_2 e^{-\varkappa x}, & x > a
\end{cases}
\qquad\quad\!
k : \ \sqrt{\frac{2mV_0 a^2}{\hbar^2} - (ka)^2} = - (ka) \cot(ka)
``

Справа записаны уравнения на волновое число ``k``, которое задаёт энергию по формуле
```math
E = \frac{\hbar^2 k^2}{2m} - V_0
```

Чтобы перейти к обозначениям из нашей задачи, нужно взять ``\hbar = a = m = 1``.

"""

# ╔═╡ b855ca85-3d2f-4521-a0d4-6dea89fc6e50
md"""
###### Проверим, что вне ямы решение убывает как ``\:e^{-\varkappa x},\:`` где ``\:\varkappa = \sqrt{2E}:``
"""

# ╔═╡ 9f1a3b08-d15b-4b2e-8f25-4b5f38581cfd
# Экспонента ли? 
begin
	offset = 0 		# Номер точки, отсчитывая от правой границы ямы
	Δ = 4 			# На сколько точек отступить вправо
	
	u[N+M+1+offset, n] * exp(-Δ*h√(-2*E[n])) / u[N+M+1+offset+Δ, n]
	# Реальное изменение u(1) / u(1+hΔ) разделить на теоретическое e^-ϰhΔ.
	# Должны получить приблизительно 1.
end

# ╔═╡ 74b094de-51bd-43c9-8c26-3b1add8166ca
md"""
###### Число решений
"""

# ╔═╡ 99b0ad2b-045c-4f1e-9b01-95ca682d82fc
md"""
Уравнения на ``k`` удобно решать графически, поскольку их левая часть задаёт дугу окружности. Можно показать, что существует

```math
\begin{array}{cl}
\left\lceil \dfrac{\sqrt{2V_0}}{\pi} \right\rceil
\qquad &\text{ чётных решений} \\ 
\left\lceil \dfrac{\sqrt{2V_0} - \pi /2 }{\pi} \right\rceil
\qquad &\text{ нечётных решений}
\end{array}
```
"""

# ╔═╡ 47513d45-f064-4165-a636-2ec4a9fe455f
# Теоретическое количество решений с отрицательной энергией
ceil(√(2V₀)/π) + ceil((√(2V₀) - π/2)/π)

# ╔═╡ 4ffcd4bf-90e1-467b-9087-6b048f210b24
# Эксперименатальное число решений
findlast(signbit, E)
# Энергии расположены по возрастанию, 
# функция signbit возвращает true, когда число отрицательное
# функция findlast возвращает индекс последнего отрицательного элемента

# ╔═╡ 7815ed79-a58d-4b71-ba12-4606cdeb1380
# Верхняя граница на энергию основного состояния
(π/2)^2/2 - V₀

# ╔═╡ 7b033426-c49e-441c-b0e5-81c77a0da8e9
md"""
##### Поведение при изменении числа элементов
"""

# ╔═╡ eadbcaee-f8c5-4d2e-90d5-18899d9316bf
md"""
Будем варьировать ``\,h = 1/M:``
"""

# ╔═╡ b27faac9-bd33-4242-99ec-257455ca9d3d
M_List = [8, 16, 32, 64, 128] 

# ╔═╡ 38920525-e37c-496b-ac4b-275dcfe31a34
n_sol = 4

# ╔═╡ e35e1b53-a2c7-4b66-b46c-0f6e7eb5de9d
EnergiesForDifferentNumberOfElements = Array{Float64}(undef, length(M_List), n_sol);

# ╔═╡ d9a6a83c-c16a-4615-bdc5-2244d44d7e6a
md"""
При каждом ``M`` решим задачу и запомним первые `n_sol = `$n_sol энергий.

Затем построим график зависимости уровней энергии от ``M``.
"""

# ╔═╡ 40522e55-533e-423a-9679-db5e0a172cf6
for i ∈ 1 : length(M_List)
	_E, = shrodingerSolve(;params..., M = M_List[i])
	EnergiesForDifferentNumberOfElements[i, :] = _E[1:n_sol]
end

# ╔═╡ fe663a3e-1ee2-4537-95a1-43d73c25a95c
plot(M_List, EnergiesForDifferentNumberOfElements,
	xlabel = "M (число элементов на единицу длины)",
	ylabel = "Уровни энергии",
	label = reshape(["E<sub>$k</sub>" for k ∈ 1 : n_sol], (1, :))
)

# ╔═╡ Cell order:
# ╠═bff05cbb-0331-47c5-be77-d8c3427eb8d6
# ╟─5f3959b0-865c-11ec-13e2-e72527971dcb
# ╟─b7667130-f054-4962-8937-8bd50b0aab4c
# ╟─d2e5fed7-1f31-4dac-8a25-1166ad1aa2d0
# ╟─00b7ca39-7a84-4a44-ac09-02e587b1aae5
# ╟─6b1e8f85-e38f-435f-b865-0e505aa16bd2
# ╟─166fd752-34c1-41b9-8032-26c86a7ba725
# ╟─af310304-aa61-49ef-a1db-90c8abd89dbf
# ╟─802591f2-c124-4663-9c30-b0ebc0c0c2ce
# ╟─85e2e1dd-c618-43d8-bedd-022ebc86ee24
# ╟─a9152c5a-30c1-4120-a8ff-523b10a8b645
# ╠═30efff26-9f0e-4ef3-8c3d-668a44e8c849
# ╠═9806ac50-5f3e-470e-a5ef-6b33b285712d
# ╠═3d9ed0c9-5f9e-4514-9593-8637c466177c
# ╠═e761b552-bdce-4abe-a6aa-f12a1302fc3e
# ╠═23a1d88c-756a-46e3-9ea4-d70a9260ebc7
# ╟─8b433c95-8b23-4c05-881a-fc53249472b6
# ╠═bcbc3be4-6414-4d9d-ad12-e28edc25844a
# ╟─d03d4e6e-4f14-41e4-a37c-a345783ea479
# ╟─b855ca85-3d2f-4521-a0d4-6dea89fc6e50
# ╠═9f1a3b08-d15b-4b2e-8f25-4b5f38581cfd
# ╟─74b094de-51bd-43c9-8c26-3b1add8166ca
# ╟─99b0ad2b-045c-4f1e-9b01-95ca682d82fc
# ╠═47513d45-f064-4165-a636-2ec4a9fe455f
# ╠═4ffcd4bf-90e1-467b-9087-6b048f210b24
# ╠═7815ed79-a58d-4b71-ba12-4606cdeb1380
# ╟─7b033426-c49e-441c-b0e5-81c77a0da8e9
# ╟─eadbcaee-f8c5-4d2e-90d5-18899d9316bf
# ╠═b27faac9-bd33-4242-99ec-257455ca9d3d
# ╠═38920525-e37c-496b-ac4b-275dcfe31a34
# ╠═e35e1b53-a2c7-4b66-b46c-0f6e7eb5de9d
# ╟─d9a6a83c-c16a-4615-bdc5-2244d44d7e6a
# ╠═40522e55-533e-423a-9679-db5e0a172cf6
# ╠═fe663a3e-1ee2-4537-95a1-43d73c25a95c
