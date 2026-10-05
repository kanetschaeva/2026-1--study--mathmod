using SciMLBase: ODEProblem, solve
using OrdinaryDiffEq: Tsit5
using Plots

# --------------------------------------------------
# Лабораторная работа №6
# Задача об эпидемии
# Вариант 2
# --------------------------------------------------

# Общая численность населения
N = 25000.0

# Начальные условия варианта
I0 = 150.0
R0 = 15.0

# Число восприимчивых
S0 = N - I0 - R0

# Коэффициенты модели
# alpha - коэффициент заболеваемости
# beta  - коэффициент выздоровления
alpha = 0.000005
beta = 0.1

# Временной интервал
t0 = 0.0
tmax = 400.0
dt = 0.1

tspan = (t0, tmax)
save_times = t0:dt:tmax

# Начальный вектор
# u[1] = S
# u[2] = I
# u[3] = R
u0 = [S0, I0, R0]

# Каталог результатов
mkpath("plots/01_epidemic")

println("="^60)
println("ЛАБОРАТОРНАЯ РАБОТА №6")
println("Задача об эпидемии")
println("Вариант 2")
println("="^60)

println()
println("Начальные условия:")
println("N = ", N)
println("S(0) = ", S0)
println("I(0) = ", I0)
println("R(0) = ", R0)

println()
println("Коэффициенты:")
println("alpha = ", alpha)
println("beta = ", beta)

# --------------------------------------------------
# Модель эпидемии
#
# Если I > I*:
#
# dS/dt = -alpha*S*I
# dI/dt =  alpha*S*I - beta*I
# dR/dt =  beta*I
#
# Если I <= I*:
#
# dS/dt = 0
# dI/dt = -beta*I
# dR/dt =  beta*I
# --------------------------------------------------

function epidemic!(du, u, p, t)

    S = u[1]
    I = u[2]
    R = u[3]

    if I > p.Istar

        # Инфицированных больше критического уровня:
        # инфекция распространяется

        du[1] = -p.alpha * S * I
        du[2] =  p.alpha * S * I - p.beta * I
        du[3] =  p.beta * I

    else

        # Инфицированных не больше критического уровня:
        # больные изолированы

        du[1] = 0.0
        du[2] = -p.beta * I
        du[3] =  p.beta * I

    end
end

# ==================================================
# СЛУЧАЙ 1
# I(0) <= I*
#
# I(0) = 150
# I* = 200
# ==================================================

Istar1 = 200.0

params1 = (
    alpha = alpha,
    beta = beta,
    Istar = Istar1
)

prob1 = ODEProblem(
    epidemic!,
    u0,
    tspan,
    params1
)

sol1 = solve(
    prob1,
    Tsit5();
    saveat = save_times,
    abstol = 1e-8,
    reltol = 1e-8
)

S1 = [u[1] for u in sol1.u]
I1 = [u[2] for u in sol1.u]
R1 = [u[3] for u in sol1.u]

println()
println("="^60)
println("СЛУЧАЙ 1: I(0) <= I*")
println("="^60)

println("I(0) = ", I0)
println("I* = ", Istar1)

println(
    "Максимальное число инфицированных: ",
    round(maximum(I1); digits=2)
)

println(
    "S в конце: ",
    round(S1[end]; digits=2)
)

println(
    "I в конце: ",
    round(I1[end]; digits=4)
)

println(
    "R в конце: ",
    round(R1[end]; digits=2)
)

# График первого случая

p1 = plot(
    sol1.t,
    S1,
    label = "S(t) - восприимчивые",
    xlabel = "Время",
    ylabel = "Число людей",
    title = "Случай 1: I(0) <= I*",
    lw = 2,
    grid = true
)

plot!(
    p1,
    sol1.t,
    I1,
    label = "I(t) - инфицированные",
    lw = 2
)

plot!(
    p1,
    sol1.t,
    R1,
    label = "R(t) - иммунитет",
    lw = 2
)

savefig(
    p1,
    "plots/01_epidemic/case1_isolation.png"
)

# ==================================================
# СЛУЧАЙ 2
# I(0) > I*
#
# I(0) = 150
# I* = 100
# ==================================================

Istar2 = 100.0

params2 = (
    alpha = alpha,
    beta = beta,
    Istar = Istar2
)

prob2 = ODEProblem(
    epidemic!,
    u0,
    tspan,
    params2
)

sol2 = solve(
    prob2,
    Tsit5();
    saveat = save_times,
    abstol = 1e-8,
    reltol = 1e-8
)

S2 = [u[1] for u in sol2.u]
I2 = [u[2] for u in sol2.u]
R2 = [u[3] for u in sol2.u]

peak_index = argmax(I2)
peak_I = I2[peak_index]
peak_t = sol2.t[peak_index]

println()
println("="^60)
println("СЛУЧАЙ 2: I(0) > I*")
println("="^60)

println("I(0) = ", I0)
println("I* = ", Istar2)

println(
    "Пик числа инфицированных: ",
    round(peak_I; digits=2)
)

println(
    "Время пика: ",
    round(peak_t; digits=2)
)

println(
    "S в конце: ",
    round(S2[end]; digits=2)
)

println(
    "I в конце: ",
    round(I2[end]; digits=4)
)

println(
    "R в конце: ",
    round(R2[end]; digits=2)
)

# График второго случая

p2 = plot(
    sol2.t,
    S2,
    label = "S(t) - восприимчивые",
    xlabel = "Время",
    ylabel = "Число людей",
    title = "Случай 2: I(0) > I*",
    lw = 2,
    grid = true
)

plot!(
    p2,
    sol2.t,
    I2,
    label = "I(t) - инфицированные",
    lw = 2
)

plot!(
    p2,
    sol2.t,
    R2,
    label = "R(t) - иммунитет",
    lw = 2
)

savefig(
    p2,
    "plots/01_epidemic/case2_epidemic.png"
)

# --------------------------------------------------
# Сравнение двух случаев
# --------------------------------------------------

p_all = plot(
    p1,
    p2;
    layout = (1, 2),
    size = (1400, 500)
)

savefig(
    p_all,
    "plots/01_epidemic/both_cases.png"
)

println()
println("="^60)
println("РАСЧЁТ ЗАВЕРШЁН")
println("="^60)

println()
println("Сохранены графики:")
println("plots/01_epidemic/case1_isolation.png")
println("plots/01_epidemic/case2_epidemic.png")
println("plots/01_epidemic/both_cases.png")
