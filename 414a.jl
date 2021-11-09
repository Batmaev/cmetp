### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 3fd4daff-9a0e-4d45-b7e5-e36f3f7e73bb
begin
	import Pkg
	Pkg.activate(mktempdir())
	
	Pkg.add("PlutoUI")
	using PlutoUI
	
	Pkg.add("Plots")
	using Plots
	plotly()
	
	Pkg.add("Roots")
	import Roots
end

# ╔═╡ 13bff350-31d3-11ec-3b13-f1bbb919f886
md"""

### Одномерная газовая динамика
### Распад произвольного разрыва
### Задание 4, №10 по книге
"""

# ╔═╡ 2444095d-48bd-4acb-b77d-02a30ca218a1
md"""
В момент времени ``t = 0`` среда в полупространстве ``x < 0`` однородна и характеризуется плотностью ``\rho_1``, давлением ``p_1`` и скоростью ``u_1``.

Аналогично в полупространстве ``x > 0``.

Требуется определить, что будет при ``t > 0``.

Задача симметрична относительно сдвига в плоскости ``yOz``.
"""

# ╔═╡ 2aa8688e-ec03-4e02-9041-5aa8661d8b13
md"""
Введём векторы

$$\vec{f} = \begin{pmatrix}
\rho \\ U \\ E 
\end{pmatrix},
\qquad
\vec{F} = \begin{pmatrix} 
U \\ U^2/\rho + p \\ (U/\rho)(E+p)
\end{pmatrix}$$

где ``\:U = \rho u\:`` -- импульс единицы объёма, \
``E = \dfrac{\rho u^2}{2} + E_{внутр}\:`` -- энергия единицы объёма, \
``E_{внутр} = \dfrac{\frac{i}{2}pV}{V} = \dfrac1{\gamma - 1}\, p``, \
``i`` -- число степеней свободы, ``\gamma`` -- показатель адиабаты.

Тогда уравнение динамики запишется в виде

$$\frac{\partial \vec f}{\partial t} + \frac{\partial \vec F}{\partial x} = 0$$

В качестве граничных условий выберем

$$\vec f|_{x=a} = \vec f(0,a) \qquad
\vec f|_{x=b} = \vec f(0,b)$$
"""

# ╔═╡ dbc58385-6d41-43c0-80eb-392898ae36de
md"""
###### Будем решать его методом Мак-Кормака.

Пусть ``\grave f`` -- текущий временной слой.

Тогда промежуточный временно́й слой:

``\displaystyle
\bar f_i = \grave f_i - \frac{\tau}{h}(\grave F_{i+1} - \grave F_i)``

И следующий временно́й слой:

`` \displaystyle
\hat f_i = \frac12 \cdot (\grave f_i + \bar f_n) - \frac{\tau}{2h} (\bar F_i - \bar F_{i-1})``


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


###### Сглаживание:

Вычислим разности плотностей в соседних узлах:

``D_i = \hat \rho_{i+1} - \hat \rho_i``

Затем вычислим векторные потоки ``Q``:

$$\begin{aligned}
Q_i^{plus} &= D_i D_{i-1} \leqslant 0 \;\;\lor\;\; D_i D_{i+1} \leqslant 0 &&? \quad \hat f_{i+1} - \hat f_i \quad : \quad 0
\\
Q_i^{minus} &= D_i D_{i-1} \leqslant 0 \;\;\lor\;\; D_{i-1} D_{i-2} \leqslant 0 &&? \quad \hat f_i - \hat f_{i-1} \quad : \quad 0
\end{aligned}$$

После этого обновляем значения на новом временном слое в тех точках, где это возможно:

$$\hat f_i \mathrel{+}= q \cdot (Q^p_i - Q^m_i), \qquad q \in [0.05, 0.2] - \text{ произвольное}$$

"""

# ╔═╡ 290ca450-cb5f-49e1-9a83-7950796f09bb
# Параметры задачи:

begin
	ρ₁ = 1 			# Плотность газа, кг/м³
	ρ₂ = 2
	
	p₁ = 50000 		# Давление, Па
	p₂ = 75000
	
	u₁ = 50 		# Скорости газов, м/c
	u₂ = -50
	
	γ = 5/2 		# Показатель адиабаты
end;

# ╔═╡ 2c64f0a4-b6fe-4c19-88cf-c337c3b8b8c0
# Параметры схемы:

begin
	a = -10 		# Левая граничная точка, м
	b = 10 			# Правая
	Nx = 1000 		# Число шагов по пространству
					# (число точек с учётом граничных будет Nx+1)
	h = (b-a) / Nx 	# Шаг по пространству, м
	
	T = 0.022 	 	# Время моделирования
	Nt = 5000
	τ = T / Nt 		
	
	σ = τ/h 		# Число Куранта
	
	q = .05 		# Параметр сглаживания
end;

# ╔═╡ 65967e0b-b1b7-4230-92d9-081a81264e01
params = (; ρ₁, ρ₂, p₁, p₂, u₁, u₂, γ, a, b, Nx, Nt, τ, q);

# ╔═╡ 16c1556d-c828-43d1-a8e3-b20c46df1b6f
md"``\tau/h`` = $σ"

# ╔═╡ 324add79-621a-46ff-907f-f50539f1110c
function diffur(;ρ₁, ρ₂, p₁, p₂, u₁, u₂, γ, a, b, Nx, Nt, τ, q)
	
	h = (b-a) / Nx
	σ = τ/h
	
	# Здесь будет храниться решение
	f = Array{Float64}(undef, Nx+1, 3, Nt+1)
	p = Array{Float64}(undef, Nx+1, Nt+1)
	
	# Первый временной слой
	ρ̀ = @view f[:, 1, 1]
	Ù = @view f[:, 2, 1]
	È = @view f[:, 3, 1]

	for (i, x) ∈ zip(1:Nx+1, a:h:b)
		if x < 0
			ρ̀[i] = ρ₁
			Ù[i] = ρ₁ * u₁
			È[i] = ρ₁*u₁^2/2 + 1/(γ-1)*p₁
		else
			ρ̀[i] = ρ₂
			Ù[i] = ρ₂ * u₂
			È[i] = ρ₂*u₂^2/2 + 1/(γ-1)*p₂
		end
	end
	
	# Будут перезаписываться на каждой итерации:
	f̄ = Array{Float64}(undef, Nx+1, 3)
	p̄ = Array{Float64}(undef, Nx+1)
	F̀ = Array{Float64}(undef, Nx+1, 3)
	F̄ = Array{Float64}(undef, Nx+1, 3)
	D = Array{Float64}(undef, Nx) 	# Разности плотностей в соседних точках
	Qᵖ= Array{Float64}(undef, 3)
	Qᵐ= Array{Float64}(undef, 3)


	
	# Основной цикл
	
	for t = 1:Nt
		
	# Текущий временной слой
		f̀ = @view f[:, :, t]
		ρ̀ = @view f[:, 1, t]
		Ù = @view f[:, 2, t]
		È = @view f[:, 3, t]
		p̀ = @view p[:, t]
		
		@. p̀ = (È - Ù^2 / 2ρ̀) * (γ-1) 
		
		@. F̀[:, 1] = Ù
		@. F̀[:, 2] = Ù^2 / ρ̀ + p̀
		@. F̀[:, 3] = Ù/ρ̀ * (È + p̀)
		

		
	# Промежуточный временной слой
		for i ∈ 2:Nx
			@. f̄[i, :] = f̀[i, :] - σ*(F̀[i+1, :] - F̀[i, :])
		end
		
		@. f̄[1, :] = f̀[1, :]
		@. f̄[Nx+1, :] = f̀[Nx+1, :]
		
		# Подсчёт F̄ на промежуточном временном слое
		ρ̄ = @view f̄[:, 1]
		Ū = @view f̄[:, 2]
		Ē = @view f̄[:, 3]
		
		@. p̄ = (Ē - Ū^2 / 2ρ̄) * (γ-1) 
		
		@. F̄[:, 1] = Ū
		@. F̄[:, 2] = Ū^2 / ρ̄ + p̄
		@. F̄[:, 3] = Ū/ρ̄ * (Ē + p̄)

		
	# Новый временной слой
		f̂ = @view f[:, :, t+1]
		
		for i ∈ 2:Nx
			@. f̂[i, :] = 1/2*(f̀[i, :] + f̄[i, :]) - σ/2*(F̄[i, :] - F̄[i-1, :])
		end
		
		@. f̂[1, :] = f̀[1, :]
		@. f̂[Nx+1, :] = f̀[Nx+1, :]
		
	# Сглаживание
		ρ̂ = @view f̂[:, 1]
		
		for i ∈ 1:Nx
			D[i] = ρ̂[i+1] - ρ̂[i]
		end
		
		for i ∈ 3:Nx-1
			Qᵖ .= D[i]D[i-1] ≤ 0 || D[i  ]D[i+1] ≤ 0 ? (@views f̂[i+1,:]-f̂[i  , :]) : 0
			Qᵐ .= D[i]D[i-1] ≤ 0 || D[i-1]D[i-2] ≤ 0 ? (@views f̂[i  ,:]-f̂[i-1, :]) : 0
			@. f̂[i, :] += q*(Qᵖ - Qᵐ)
		end

	end
	
	return f, p
	
end;

# ╔═╡ 01f33e2e-e9b7-490f-b91d-6747a5f9811f
4 * (Nt * Nx) * 8 /1024/1024 # Megabytes

# ╔═╡ 4e6b11d2-2f71-4a53-9b88-3108d8ed2996
f, p = diffur(; params...);

# ╔═╡ b4730cb1-b459-4ba0-b26c-df83b6430090
@bind t Slider(1:Nt)

# ╔═╡ 8c2d07d4-adfd-4bb4-a759-37cc5301e29b
md"""
##### Аналитическое решение

Решение определяется скоростью сближения газов: ``\; u_1 - u_2``

Есть три критические скорости сближения:
"""

# ╔═╡ e1e8e63f-5f58-4f53-9ea5-e23e801c84b4
Uₛₕ = (p₂ - p₁) / √(ρ₁*((γ+1)p₂ + (γ-1)p₁)/2)

# ╔═╡ ed0c48e2-c8e1-4b1f-b9f7-d632acf05764
begin
	c₁ = √(γ*p₁/ρ₁)
	c₂ = √(γ*p₂/ρ₂)
end;

# ╔═╡ a7e500bc-bd04-4812-81b5-c80501547771
Uₒᵤₜ = -2c₂/(γ-1) * (1 - (p₁/p₂)^((γ-1)/2γ))

# ╔═╡ 44b727a5-ee7b-4f72-8e87-1375a818e7e7
Uᵥₐₖ = -2(c₁ + c₂) / (γ-1)

# ╔═╡ 21a026dc-94cb-4ede-861b-4ed96212aaee
md"""

|       ``\quad u_1 - u_2``        | Конфигурация решения |
|:--------------------------|:---------------------|
|``\in (-\infty,\, U_{vak})``| Вправо и влево распространяются волны разрежения, в центре область вакуума|
|``\in (U_{vak},\, U_{out})``| Вправо и влево распространяются волны разрежения, в центре область пониженного давления|
|``\in (U_{out},\, U_{sh})`` | Налево распространяется ударная волна, направо -- волна разрежения|
|``\in (U_{sh}, \, +\infty)``| Вправо и влево распространяются ударные волны|

Вычислим по готовым формулам аналитическое решение в случае двух ударных волн.

"""

# ╔═╡ 38dd3b3e-fe52-47c6-936e-df05763132f5
begin
	a₁(p) = √(ρ₁*((γ+1)p + (γ-1)p₁)/2)
	a₂(p) = √(ρ₂*((γ+1)p + (γ-1)p₂)/2)

	φ(p) = (a₂(p)*p₁ + a₁(p)*p₂ + a₁(p)*a₂(p)*(u₁-u₂)) / (a₁(p)+a₂(p))
end

# ╔═╡ 94ad644d-813c-4e95-8e95-60950e6b4eb4
# Давление посередине
P = Roots.find_zero(p -> φ(p) - p, (p₁+p₂)/2)

# ╔═╡ 111b1da9-1518-4bd4-8bb6-b34935b74877
# Скорость границы между левой и правой волной
u = (a₁(P)u₁ + a₂(P)u₂ + p₁ - p₂) / (a₁(P) + a₂(P))

# ╔═╡ 21620d93-6575-44f6-a78b-ee32477e0de9
# Плотность в левой волне
R₁ = ρ₁ * ((γ+1)P + (γ-1)p₁) / ((γ-1)P + (γ+1)p₁)

# ╔═╡ 130763bf-4ab8-4e71-8e08-83614002dbc9
# Плотность в правой волне
R₂ = ρ₂ * ((γ+1)P + (γ-1)p₂) / ((γ-1)P + (γ+1)p₂)

# ╔═╡ 6e2271a4-e403-446a-a68d-7b0a8bcdeaee
# Скорость левой волны
D₁ = u₁ - a₁(P) / ρ₁

# ╔═╡ 5942ddc9-bb88-4a4f-a69f-b1576cfacdbd
# Скорость правой волны
D₂ = u₂ + a₂(P) / ρ₂

# ╔═╡ f9d704d3-a106-4a2f-a573-1ca0ce12db09
anal = Array{Float64}(undef, Nx+1, 2, Nt+1);

# ╔═╡ e99f93d5-ccb4-49b6-b36f-44e33c92fef4
begin
	plot(a:h:b, [f[:, 1, t], p[:, t]],
		layout = (2,1),
		label = ["Численное" :none],
		legend = :bottomright,
		ylabel = ["ρ, кг/м³" "P, Па"],
		title = ["t = $(round((t-1)τ, digits = 5)) с" ""],
		xlabel = ["" "x, м"],
	)
	try
		plot!(a:h:b, anal[:, :, t], label = ["Аналитическое" :none],)
	catch
	end
	plot!()
end

# ╔═╡ 80146061-6900-450a-b94d-6fdefc214cd1
for (i_t, t) ∈ zip(1:Nt, 0:τ:T)
	for (i_x, x) ∈ zip(1:Nx+1, a:h:b)
		anal[i_x, :, i_t] = 
			if x < D₁*t
				[ρ₁, p₁]
			elseif x < u*t
				[R₁, P]
			elseif x < D₂*t
				[R₂, P]
			else
				[ρ₂, p₂]
			end
	end
end

# ╔═╡ 08b0a22e-3c65-475e-9eb0-9a93cf69d489
md"""
##### Численное решение при других конфигурациях
"""

# ╔═╡ 779fb175-6721-424d-9b80-81b436237d83
md"""
Будем верхним индексом обозначать номер конфигурации.

###### 1). Влево распространяется ударная волна, вправо -- волна разрежения.
"""

# ╔═╡ 80902f6e-b8ba-484f-87a6-fc65ff8b54df
begin
	Δu¹ = (Uₒᵤₜ + Uₛₕ) / 2  	# Возьмём разность скоростей из середины допустимого промежутка
	u₁¹ = Δu¹ / 3
	u₂¹ = -2/3 * Δu¹
end;

# ╔═╡ 94faa92a-ba5e-42fb-8e57-7fd38e132266
(u₁¹, u₂¹)

# ╔═╡ 05f8775b-6167-4763-b043-14a2a1837f7c
f¹, p¹ = diffur(; params..., u₁=u₁¹, u₂=u₂¹);

# ╔═╡ 0b5b18b0-7d37-40da-859e-f888242f8a4e
@bind t¹ Slider(1:Nt)

# ╔═╡ 357b2085-755e-4254-9e27-b5c30fa1ac8d
plot(a:h:b, [f¹[:, 1, t¹], p¹[:, t¹]],
	layout = (2,1),
	label = false,
	ylabel = ["ρ, кг/м³" "P, Па"],
	title = ["t = $(round((t¹-1)τ, digits = 5)) с" ""],
	xlabel = ["" "x, м"],
)

# ╔═╡ 98f98623-e69e-4ab8-8698-ca08d5224028
md"###### 2). Две волны разрежения"

# ╔═╡ 226946f3-8047-4faf-9e7a-348066d28dac
begin
	Δu² = (Uᵥₐₖ + Uₒᵤₜ) / 2  	# Возьмём разность скоростей из середины допустимого промежутка
	u₁² =  Δu² / 2
	u₂² = -Δu² / 2
end;

# ╔═╡ 1e5786a2-f579-4c50-abf6-8d73f23faa15
(u₁², u₂²)

# ╔═╡ 1ae2ff8d-349b-4c28-87c1-6ece78fe1273
begin
	Nt² = 2Nt÷3
	Nx² = Nx
	h² = (b-a)/Nx²
	τ² = τ
end

# ╔═╡ 376555e5-c0db-426d-bc54-612236ece9fb
4 * (Nt² * Nx²) * 8 /1024/1024 # Megabytes

# ╔═╡ 4f13698e-037c-4b47-ae7f-d1146e1c6e43
f², p² = diffur(; params..., u₁ = u₁², u₂ = u₂², Nt = Nt², τ = τ², Nx = Nx²);

# ╔═╡ 2bbe14a2-5897-4f50-8863-f8c1de15eb63
@bind t² Slider(1:Nt²)

# ╔═╡ 982006d5-8047-48f4-87ec-0e5ca3bc0491
plot(a:h²:b, [f²[:, 1, t²], p²[:, t²]],
	layout = (2,1),
	label = false,
	ylabel = ["ρ, кг/м³" "P, Па"],
	title = ["t = $(round((t²-1)τ², digits = 5)) с" ""],
	xlabel = ["" "x, м"],
)

# ╔═╡ cc84ce4e-f6a7-4360-8128-fe23142b813d
md"""###### 3). Вакуум в центре"""

# ╔═╡ 44b58e5a-c014-4b8d-bd3f-8d3ac34f4495
begin
	u₁³ = 3/4 * Uᵥₐₖ
	u₂³ = -3/4 * Uᵥₐₖ
end;

# ╔═╡ b8300ba0-32dc-4aeb-a0fa-c164b1f8e931
(u₁³, u₂³)

# ╔═╡ 7499221a-9a80-408c-881e-82663d113045
τ

# ╔═╡ 1d6913da-4ae3-4ba9-8987-08c8f2bc0034
begin
	τ³ = τ/10
	Nt³ = 2Nt
	Nx³ = 2Nx
	h³ = (b-a) / Nx³
end;

# ╔═╡ 38e7fa40-f50d-4d8a-98a8-c35af13e8c0c
4 * (Nx³+1)*Nt³ * 8 / 1024/1024 # Megabytes

# ╔═╡ ef394459-2ba3-4ba3-ac3c-79adc2aea609
f³, p³ = diffur(; params..., u₁=u₁³, u₂=u₂³, τ = τ³, Nt = Nt³, Nx = Nx³);

# ╔═╡ 0c06f6cb-6b2b-43e1-9f60-bb9d01ae86a7
@bind t³ Slider(1:Nt³)

# ╔═╡ 08a87d83-da5c-42da-a307-114a67079fee
plot(a:h³:b, [f³[:, 1, t³], p³[:, t³]],
	layout = (2,1),
	label = false,
	ylabel = ["ρ, кг/м³" "P, Па"],
	title = ["t = $(round((t³-1)τ³, digits = 5)) с" ""],
	xlabel = ["" "x, м"],
)

# ╔═╡ Cell order:
# ╠═3fd4daff-9a0e-4d45-b7e5-e36f3f7e73bb
# ╟─13bff350-31d3-11ec-3b13-f1bbb919f886
# ╟─2444095d-48bd-4acb-b77d-02a30ca218a1
# ╟─2aa8688e-ec03-4e02-9041-5aa8661d8b13
# ╟─dbc58385-6d41-43c0-80eb-392898ae36de
# ╠═290ca450-cb5f-49e1-9a83-7950796f09bb
# ╠═2c64f0a4-b6fe-4c19-88cf-c337c3b8b8c0
# ╠═65967e0b-b1b7-4230-92d9-081a81264e01
# ╟─16c1556d-c828-43d1-a8e3-b20c46df1b6f
# ╠═324add79-621a-46ff-907f-f50539f1110c
# ╠═01f33e2e-e9b7-490f-b91d-6747a5f9811f
# ╠═4e6b11d2-2f71-4a53-9b88-3108d8ed2996
# ╟─b4730cb1-b459-4ba0-b26c-df83b6430090
# ╠═e99f93d5-ccb4-49b6-b36f-44e33c92fef4
# ╟─8c2d07d4-adfd-4bb4-a759-37cc5301e29b
# ╟─e1e8e63f-5f58-4f53-9ea5-e23e801c84b4
# ╟─a7e500bc-bd04-4812-81b5-c80501547771
# ╟─44b727a5-ee7b-4f72-8e87-1375a818e7e7
# ╟─ed0c48e2-c8e1-4b1f-b9f7-d632acf05764
# ╟─21a026dc-94cb-4ede-861b-4ed96212aaee
# ╠═94ad644d-813c-4e95-8e95-60950e6b4eb4
# ╠═38dd3b3e-fe52-47c6-936e-df05763132f5
# ╠═111b1da9-1518-4bd4-8bb6-b34935b74877
# ╠═21620d93-6575-44f6-a78b-ee32477e0de9
# ╠═130763bf-4ab8-4e71-8e08-83614002dbc9
# ╠═6e2271a4-e403-446a-a68d-7b0a8bcdeaee
# ╠═5942ddc9-bb88-4a4f-a69f-b1576cfacdbd
# ╠═f9d704d3-a106-4a2f-a573-1ca0ce12db09
# ╠═80146061-6900-450a-b94d-6fdefc214cd1
# ╟─08b0a22e-3c65-475e-9eb0-9a93cf69d489
# ╟─779fb175-6721-424d-9b80-81b436237d83
# ╠═80902f6e-b8ba-484f-87a6-fc65ff8b54df
# ╠═94faa92a-ba5e-42fb-8e57-7fd38e132266
# ╠═05f8775b-6167-4763-b043-14a2a1837f7c
# ╟─0b5b18b0-7d37-40da-859e-f888242f8a4e
# ╟─357b2085-755e-4254-9e27-b5c30fa1ac8d
# ╟─98f98623-e69e-4ab8-8698-ca08d5224028
# ╠═226946f3-8047-4faf-9e7a-348066d28dac
# ╠═1e5786a2-f579-4c50-abf6-8d73f23faa15
# ╠═1ae2ff8d-349b-4c28-87c1-6ece78fe1273
# ╠═376555e5-c0db-426d-bc54-612236ece9fb
# ╠═4f13698e-037c-4b47-ae7f-d1146e1c6e43
# ╟─2bbe14a2-5897-4f50-8863-f8c1de15eb63
# ╟─982006d5-8047-48f4-87ec-0e5ca3bc0491
# ╟─cc84ce4e-f6a7-4360-8128-fe23142b813d
# ╠═44b58e5a-c014-4b8d-bd3f-8d3ac34f4495
# ╠═b8300ba0-32dc-4aeb-a0fa-c164b1f8e931
# ╠═7499221a-9a80-408c-881e-82663d113045
# ╠═1d6913da-4ae3-4ba9-8987-08c8f2bc0034
# ╠═38e7fa40-f50d-4d8a-98a8-c35af13e8c0c
# ╠═ef394459-2ba3-4ba3-ac3c-79adc2aea609
# ╟─0c06f6cb-6b2b-43e1-9f60-bb9d01ae86a7
# ╠═08a87d83-da5c-42da-a307-114a67079fee
