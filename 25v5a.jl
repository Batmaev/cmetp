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

# ╔═╡ 59db45bf-c9f7-4989-b92a-f0862bd7a474
begin
	import Pkg
	
	Pkg.add("PlutoUI")
	using PlutoUI
	
	Pkg.add("Plots")
	using Plots
	plotly()
	
	Pkg.add("DataFrames")
	using DataFrames
	
	Pkg.add("GLM")
	import GLM
end

# ╔═╡ efbb1b00-1e4b-11ec-32e0-3b41bf2b8b36
md"""

### Задание 2, № 5 по книге
###### Постановка задачи:

Уравнение:

``
\displaystyle
\frac{\partial u}{\partial t} + \frac{\partial F}{\partial x}
= \nu \frac{\partial^2 u}{\partial x^2},
\qquad F = \frac{u^2}{2}
``

Начальное условие: 

a) ``\; \displaystyle 
u|_{t=0} = \operatorname{HeavisideTheta}(-x) 
``

b) ``\; \displaystyle
u|_{t=0} = \exp\!\left(
	-\frac{(x-x_0)^2}{d^2}
\right)
``

Краевые условия:

`` u(t, a) \equiv u(0, a) \qquad u(t,b) \equiv u(0,b) ``

Числа Куранта:

`` \displaystyle \mu = \frac{\nu \tau}{h^2}, 
\qquad r = \frac{\tau}{h}``


"""

# ╔═╡ 0656c12b-d327-47dc-9107-2ed74f1f2623
md"""
###### Зададим параметры и начальные функции:
"""

# ╔═╡ 9a2596d1-493a-45f0-b8bf-7f379d667d58
begin
	a = -1           # Границы промежутка по пространству
	b = 1
	
	T = 1.5          # Граница промежутка по времени

	Nx = 1000        # Точки нумеруются от 1 до nX+1
	Nt = 15000
	
	τ = T / Nt
	
	h = (b - a) / Nx
	
	
	ν = 0.001 		 # Коэффициент в уравнении
	
	r = τ / h        # Гиперболическое число Куранта
	
	μ = ν * τ / h^2  # Параболическое число Куранта
	
	
	params = (a=a, b, Nx, Nt, τ, h, r)
	# ν в функцию передаём отдельно,
	# потому что его нужно будет варьировать
end;

# ╔═╡ d7f30bc9-a783-4b0e-ad1e-d60b166f3f25
u0a(x) = x > 0 ? 0. : 1.

# ╔═╡ 2ad832ed-d147-4fc9-a404-a0c4b4e13871
begin
	x0 = 0.
	d = 0.1
	
	u0b(x) = exp(-(x-x0)^2 / d^2)
end

# ╔═╡ 8730a761-5ed2-49ef-85e3-738d9bdbff40
md"""

###### Используем схему Мак-Кормака:

``\displaystyle
r = \frac{\tau}{h},
\qquad
\mu = \frac{\nu \tau}{h^2},
\qquad
F = \frac{u^2}{2}
``

1). Промежуточный временной слой:

``\qquad 
\bar u_j = u_j - r(F_{j+1} - F_j) + \mu (u_{j+1} - 2 u_j + u_{j-1})
``

2). Следующий временной слой:

``\qquad
\hat u_j = \frac{1}{2}(u_j + \bar u_j) - \frac{r}{2}(\bar F_j - \bar F_{j-1}) + \frac{\mu}{2} (\bar u_{j+1} - 2 \bar u_j + \bar u_{j+1}) ``

"""

# ╔═╡ 974b08df-41e0-4917-af66-095c1c6a4032
md"""
###### Устойчивость:

Чтобы схема Мак-Кормака была устойчивой, необходимо одновременное выполненение двух условий:

``\displaystyle
\mu \leqslant \frac{1}{2}
\qquad\quad
r \cdot \operatorname{max} |u| \leqslant 1
``
"""

# ╔═╡ bb90c7b3-3bbd-42d3-bc05-4545c5cfa8dc
md"""
У нас:

``\mu = `` $μ

``\operatorname{max}|u| \overset{?}{=} 1 ``

``r = `` $r

"""

# ╔═╡ 204c1031-52e5-4a7d-97e1-9e2054c425dd
md"""
###### Точки возле границ:

Если обозначить настоящие узлы как `*`, а промежуточные -- как `⋅`, то шаблон можно изобразить так:

```
        *
        |
    ⋅ - ⋅ - ⋅
    |   |   |
* - * - * - * - * 
```

Кажется, что схема Мак-Кормака не позволяет определить первые две и последние две точки. Однако мы можем при помощи краевых условий доопределить полностью значения сначала на промежуточном временном слое, а затем и на настоящем временном слое.  
"""

# ╔═╡ 64311ba3-9ce5-49ca-8667-2b40459ad384
md"""
###### Код:
"""

# ╔═╡ 73e25696-0c66-4b97-b156-9a6b5e93338f
begin
	xlist = a : h : b
	tlist = 0 : τ : T
end;

# ╔═╡ c63ad377-2157-4ea8-8121-e4414bff4a6a
#begin
function calcU(u0, ν, params)
	
	a, b, Nx, Nt, τ, h, r = params
	# Нужно использовать параметры, а не глобальные переменные, 
	# потому что иначе функция долго работает.
	# 
	# Эту функцию будут запускать снова и снова, варьируя ν,
	# поэтому параметр ν передаётся отдельно от остальных
	
	μ = ν*τ/h^2
	
	
	# Здесь будет храниться решение
	u = Array{Float64}(undef, Nx+1, Nt+1)
	
	# Первый временной слой
	u[:, 1] = u0.(a:h:b)
	
	# Будут перезаписываться на каждой итерации:
	ū = Array{Float64}(undef, Nx+1)
	F = Array{Float64}(undef, Nx+1)
	F̄ = Array{Float64}(undef, Nx+1)


	
	# Основной цикл
	
	for n = 1 : Nt
		
		# Текущий временной слой
		u_n = @view u[:,n]
		@. F = u_n^2 / 2

		
		# Промежуточный временной слой
		for j = 2 : Nx
			ū[j] = u_n[j] - r*(F[j+1] - F[j]) + μ*(u_n[j+1] - 2u_n[j] + u_n[j-1])
		end
		
		ū[1] = u0(a)
		ū[Nx+1] = u0(b)
		
		@. F̄ = ū^2 / 2


		# Новый временной слой
		û = @view u[:, n+1]
		
		for j = 2 : Nx
			û[j] = ( 1/2*(u_n[j] + ū[j]) 
				   - r/2*(F̄[j] - F̄[j-1]) 
				   + μ/2*(ū[j+1] - 2ū[j] + ū[j-1]) )
		end
		û[1] = u0(a)
		û[Nx+1] = u0(b)

	end
	
	return u
	
end;

# ╔═╡ e5dd4590-73fb-4eb0-a807-24e7b505647c
u = calcU(u0a, ν, params);

# ╔═╡ e35eb3f0-6aa8-41dc-89db-3710df7e56ec
@bind t Slider(1 : Nt)

# ╔═╡ 44e78562-f35d-41e3-9fa9-b9fd0a16b800
plot(
	xlist, u[:, t],
	ylims=(0, 1.1),
	title = "t = $(round((t*τ), digits=3))",
	label="",
)

# ╔═╡ 4c930687-f0a7-4269-a035-ebbed890898e
md"""
###### Движение фронта:

В случае функции-ступеньки определим скорость фронта. 

Для этого построим график зависимости положения средней точки фронта от времени.

В качестве средней точки возьмём первый встретившийся узел сетки, в котором ``u < 0.5``.
"""

# ╔═╡ 3ade1dd3-81c2-4be2-8f19-2ee4734d3d34
begin
	u_ = calcU(u0a, ν, params)
	
	front_positions = Array{Float64}(undef, Nt+1)
	
	for t ∈ 1 : Nt+1
		
		m = findfirst(y -> y < 0.5, u_[:,t])
		x = a + h*(m-1) # Положение середины фронта
		
		front_positions[t] = x
	end	
end

# ╔═╡ 649b144e-1be2-43d8-833b-f71b31981418
plot(
	tlist, front_positions,
	title = "Движение середины фронта",
	xlabel="t",
	ylabel="x",
	label="",
)

# ╔═╡ c16e3df1-cd6e-4903-8a04-b64e4fa7a404
md"""
Это прямая, то есть фронт движется с постоянной скоростью, причём
``\displaystyle \;
v = \frac{\Delta x}{\Delta t} ``
"""

# ╔═╡ 340a0a5a-ff41-4c50-bdb2-ec056229eebe
velocity = 
	(front_positions[Nt+1] - front_positions[1]) / T

# ╔═╡ 9dcf774f-3870-4f5a-b443-fad5f5ee9ee2
md"""
###### Ширина фронта:

Определим границы фронта как точки, где значение функции равно 0.9 и 0.1.

Построим зависимость ширины фронта от коэффициента вязкости ``\nu``.

Ширину будем смотреть на последнем временном слое.

"""

# ╔═╡ bdd834f9-154d-4404-9a9a-1d414d5a7a9c
begin
	lower_bound = 0.1
	upper_bound = 0.9
end;

# ╔═╡ 45784374-7c8d-466d-92dc-3ab255d2699e
νlist = 0.001 : 0.001 : 0.02;

# ╔═╡ 8f5d5c07-b3cf-401e-983a-a972744ccf3b
begin
	front_widths = Array{Float64}(undef, length(νlist))

	for (i, ν) ∈ enumerate(νlist)
		
		u = calcU(u0a, ν, params)
		
		front_start = findfirst(y -> y < upper_bound, 
								u[:,end])
		front_end   = findfirst(y -> y < lower_bound, 
								u[:,end])
		
		front_widths[i] = h * (front_end - front_start)
	
	end
end

# ╔═╡ c4d9ba53-f111-4cfd-a671-62188c0b756e
md"""
В первом приближении зависимость ширины ``w`` от вязкости ``\mu`` можно описать прямой: 

``\displaystyle w \approx 8.8 \cdot \nu ``

Вот как были вычислены коэффициенты:
"""

# ╔═╡ 05ad13c7-0dab-40e7-b815-27f8806ea92a
frontWidths_DataFrame = DataFrame(ν = νlist, w = front_widths);

# ╔═╡ 7807beaf-1500-4930-89fc-f7a559f40966
fitted_model = GLM.lm(@GLM.formula(w ~ ν), frontWidths_DataFrame)

# ╔═╡ fce086e7-66d7-49d8-a750-62557935a40a
fitted_vals = GLM.predict(fitted_model);

# ╔═╡ 28de32c2-c0c5-4ccc-bae0-ec32fa9a3c8e
begin
	scatter(
		νlist, front_widths,
		title = "Ширина фронта",
		xlabel="Вязкость",
		label="«экспериментальные» данные",
		markerstrokewidth = 0,
		markercolor = "black",
		legend = :outerright,
	)
	plot!(νlist, fitted_vals,
		label="аппроксимация",
		color="slateblue1",
	)	
end

# ╔═╡ Cell order:
# ╠═59db45bf-c9f7-4989-b92a-f0862bd7a474
# ╟─efbb1b00-1e4b-11ec-32e0-3b41bf2b8b36
# ╟─0656c12b-d327-47dc-9107-2ed74f1f2623
# ╠═9a2596d1-493a-45f0-b8bf-7f379d667d58
# ╠═d7f30bc9-a783-4b0e-ad1e-d60b166f3f25
# ╠═2ad832ed-d147-4fc9-a404-a0c4b4e13871
# ╟─8730a761-5ed2-49ef-85e3-738d9bdbff40
# ╟─974b08df-41e0-4917-af66-095c1c6a4032
# ╟─bb90c7b3-3bbd-42d3-bc05-4545c5cfa8dc
# ╟─204c1031-52e5-4a7d-97e1-9e2054c425dd
# ╟─64311ba3-9ce5-49ca-8667-2b40459ad384
# ╠═73e25696-0c66-4b97-b156-9a6b5e93338f
# ╠═c63ad377-2157-4ea8-8121-e4414bff4a6a
# ╠═e5dd4590-73fb-4eb0-a807-24e7b505647c
# ╟─e35eb3f0-6aa8-41dc-89db-3710df7e56ec
# ╠═44e78562-f35d-41e3-9fa9-b9fd0a16b800
# ╟─4c930687-f0a7-4269-a035-ebbed890898e
# ╠═3ade1dd3-81c2-4be2-8f19-2ee4734d3d34
# ╠═649b144e-1be2-43d8-833b-f71b31981418
# ╟─c16e3df1-cd6e-4903-8a04-b64e4fa7a404
# ╠═340a0a5a-ff41-4c50-bdb2-ec056229eebe
# ╟─9dcf774f-3870-4f5a-b443-fad5f5ee9ee2
# ╠═bdd834f9-154d-4404-9a9a-1d414d5a7a9c
# ╠═45784374-7c8d-466d-92dc-3ab255d2699e
# ╠═8f5d5c07-b3cf-401e-983a-a972744ccf3b
# ╟─28de32c2-c0c5-4ccc-bae0-ec32fa9a3c8e
# ╠═c4d9ba53-f111-4cfd-a671-62188c0b756e
# ╠═05ad13c7-0dab-40e7-b815-27f8806ea92a
# ╠═7807beaf-1500-4930-89fc-f7a559f40966
# ╠═fce086e7-66d7-49d8-a750-62557935a40a
