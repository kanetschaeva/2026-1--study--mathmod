using SciMLBase: ODEProblem, solve, ContinuousCallback, terminate!
using OrdinaryDiffEq: Tsit5
using Plots

# --------------------------------------------------
# Лабораторная работа №3
# Вариант 2
# --------------------------------------------------

# Начальные численности армий
x0 = 25000.0
y0 = 13000.0

# Временной интервал
t0 = 0.0
tmax = 5.0
tspan = (t0, tmax)

# Каталог для графиков
mkpath("plots/01_lanchester")

println("="^60)
println("ЛАБОРАТОРНАЯ РАБОТА №3")
println("Вариант 2")
println("="^60)

println()
println("Начальные условия:")
println("Армия X: ", x0)
println("Армия Y: ", y0)

# --------------------------------------------------
# Остановка моделирования
# --------------------------------------------------
#
# Если численность одной из сторон становится равна нулю,
# эта сторона считается проигравшей.
#

function stop_condition(u, t, integrator)
    return min(u[1], u[2])
end

function stop_affect!(integrator)
    terminate!(integrator)
end

stop_callback = ContinuousCallback(
    stop_condition,
    stop_affect!
)

# ==================================================
# МОДЕЛЬ 1
# Боевые действия между регулярными войсками
#
# dx/dt = -0.41x - 0.83y + sin(t + 3)
# dy/dt = -0.29x - 0.63y + cos(t + 3)
# ==================================================

function regular!(du, u, p, t)
    x = u[1]
    y = u[2]

    du[1] = -0.41 * x - 0.83 * y + sin(t + 3)
    du[2] = -0.29 * x - 0.63 * y + cos(t + 3)
end

u0 = [x0, y0]

prob1 = ODEProblem(
    regular!,
    u0,
    tspan,
    nothing
)

sol1 = solve(
    prob1,
    Tsit5();
    callback = stop_callback,
    saveat = 0.01,
    abstol = 1e-8,
    reltol = 1e-8
)

x1 = [max(u[1], 0.0) for u in sol1.u]
y1 = [max(u[2], 0.0) for u in sol1.u]

println()
println("="^60)
println("МОДЕЛЬ 1: РЕГУЛЯРНЫЕ ВОЙСКА")
println("="^60)

println(
    "Время окончания моделирования: ",
    round(sol1.t[end]; digits=4)
)

println(
    "Армия X: ",
    round(x1[end]; digits=2)
)

println(
    "Армия Y: ",
    round(y1[end]; digits=2)
)

if x1[end] <= 1e-3 && y1[end] > 1e-3
    println("Победитель: армия Y")
elseif y1[end] <= 1e-3 && x1[end] > 1e-3
    println("Победитель: армия X")
else
    println("Победитель не определён")
end

# График первой модели

p1 = plot(
    sol1.t,
    x1,
    label = "Армия X",
    xlabel = "Время",
    ylabel = "Численность армии",
    title = "Модель 1 - регулярные войска",
    lw = 3,
    grid = true
)

plot!(
    p1,
    sol1.t,
    y1,
    label = "Армия Y",
    lw = 3
)

savefig(
    p1,
    "plots/01_lanchester/regular_forces.png"
)

# ==================================================
# МОДЕЛЬ 2
# Регулярные войска и партизанские отряды
#
# dx/dt = -0.33x - 0.88y + sin(t)
# dy/dt = -0.44xy - 0.77y + cos(3t)
# ==================================================

function mixed!(du, u, p, t)
    x = u[1]
    y = u[2]

    du[1] = -0.33 * x - 0.88 * y + sin(t)
    du[2] = -0.44 * x * y - 0.77 * y + cos(3 * t)
end

prob2 = ODEProblem(
    mixed!,
    u0,
    tspan,
    nothing
)

sol2 = solve(
    prob2,
    Tsit5();
    callback = stop_callback,
    saveat = 0.001,
    abstol = 1e-8,
    reltol = 1e-8
)

x2 = [max(u[1], 0.0) for u in sol2.u]
y2 = [max(u[2], 0.0) for u in sol2.u]

println()
println("="^60)
println("МОДЕЛЬ 2: РЕГУЛЯРНЫЕ ВОЙСКА И ПАРТИЗАНЫ")
println("="^60)

println(
    "Время окончания моделирования: ",
    round(sol2.t[end]; digits=4)
)

println(
    "Армия X: ",
    round(x2[end]; digits=2)
)

println(
    "Армия Y: ",
    round(y2[end]; digits=2)
)

if x2[end] <= 1e-3 && y2[end] > 1e-3
    println("Победитель: армия Y")
elseif y2[end] <= 1e-3 && x2[end] > 1e-3
    println("Победитель: армия X")
else
    println("Победитель не определён")
end

# График второй модели

p2 = plot(
    sol2.t,
    x2,
    label = "Армия X",
    xlabel = "Время",
    ylabel = "Численность армии",
    title = "Модель 2 - регулярные войска и партизаны",
    lw = 3,
    grid = true
)

plot!(
    p2,
    sol2.t,
    y2,
    label = "Армия Y",
    lw = 3
)

savefig(
    p2,
    "plots/01_lanchester/mixed_forces.png"
)

# --------------------------------------------------
# Общий рисунок
# --------------------------------------------------

p_all = plot(
    p1,
    p2;
    layout = (1, 2),
    size = (1200, 500)
)

savefig(
    p_all,
    "plots/01_lanchester/both_models.png"
)

println()
println("="^60)
println("РАСЧЁТ ЗАВЕРШЁН")
println("="^60)

println()
println("Сохранены графики:")
println("plots/01_lanchester/regular_forces.png")
println("plots/01_lanchester/mixed_forces.png")
println("plots/01_lanchester/both_models.png")
