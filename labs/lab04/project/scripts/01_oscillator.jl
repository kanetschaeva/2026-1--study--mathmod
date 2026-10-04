using SciMLBase: ODEProblem, solve
using OrdinaryDiffEq: Tsit5
using Plots

# --------------------------------------------------
# Лабораторная работа №4
# Вариант 2
# Гармонический осциллятор
# --------------------------------------------------

# Начальные условия
x0 = 1.0
y0 = 1.0

# Интервал и шаг из задания
t0 = 0.0
tmax = 40.0
dt = 0.05

tspan = (t0, tmax)
save_times = t0:dt:tmax

# Вектор начальных условий:
# u[1] = x
# u[2] = x' = y
u0 = [x0, y0]

# Каталог для графиков
mkpath("plots/01_oscillator")

println("="^60)
println("ЛАБОРАТОРНАЯ РАБОТА №4")
println("Вариант 2")
println("="^60)

println()
println("Начальные условия:")
println("x(0) = ", x0)
println("y(0) = ", y0)
println("Интервал: [", t0, "; ", tmax, "]")
println("Шаг: ", dt)

# ==================================================
# СЛУЧАЙ 1
# Без затухания и без внешней силы
#
# x'' + 3x = 0
#
# Переходим к системе:
# x' = y
# y' = -3x
# ==================================================

function oscillator1!(du, u, p, t)
    x = u[1]
    y = u[2]

    du[1] = y
    du[2] = -3.0 * x
end

prob1 = ODEProblem(
    oscillator1!,
    u0,
    tspan
)

sol1 = solve(
    prob1,
    Tsit5();
    saveat = save_times
)

x1 = [u[1] for u in sol1.u]
y1 = [u[2] for u in sol1.u]

println()
println("="^60)
println("СЛУЧАЙ 1")
println("Без затухания и без внешней силы")
println("x'' + 3x = 0")
println("="^60)

println("Минимальное x: ", round(minimum(x1); digits=4))
println("Максимальное x: ", round(maximum(x1); digits=4))

# Решение x(t)

p1_solution = plot(
    sol1.t,
    x1,
    xlabel = "t",
    ylabel = "x(t)",
    title = "Случай 1 - решение x(t)",
    label = "x(t)",
    lw = 2,
    grid = true
)

savefig(
    p1_solution,
    "plots/01_oscillator/case1_solution.png"
)

# Фазовый портрет: x' от x

p1_phase = plot(
    x1,
    y1,
    xlabel = "x",
    ylabel = "x'",
    title = "Случай 1 - фазовый портрет",
    label = false,
    lw = 2,
    grid = true
)

savefig(
    p1_phase,
    "plots/01_oscillator/case1_phase.png"
)

p1_combined = plot(
    p1_solution,
    p1_phase;
    layout = (1, 2),
    size = (1200, 450)
)

savefig(
    p1_combined,
    "plots/01_oscillator/case1_combined.png"
)

# ==================================================
# СЛУЧАЙ 2
# С затуханием, без внешней силы
#
# x'' + x' + 4x = 0
#
# Система:
# x' = y
# y' = -y - 4x
# ==================================================

function oscillator2!(du, u, p, t)
    x = u[1]
    y = u[2]

    du[1] = y
    du[2] = -y - 4.0 * x
end

prob2 = ODEProblem(
    oscillator2!,
    u0,
    tspan
)

sol2 = solve(
    prob2,
    Tsit5();
    saveat = save_times
)

x2 = [u[1] for u in sol2.u]
y2 = [u[2] for u in sol2.u]

println()
println("="^60)
println("СЛУЧАЙ 2")
println("С затуханием, без внешней силы")
println("x'' + x' + 4x = 0")
println("="^60)

println("x в конце интервала: ", round(x2[end]; digits=6))
println("x' в конце интервала: ", round(y2[end]; digits=6))

p2_solution = plot(
    sol2.t,
    x2,
    xlabel = "t",
    ylabel = "x(t)",
    title = "Случай 2 - решение x(t)",
    label = "x(t)",
    lw = 2,
    grid = true
)

savefig(
    p2_solution,
    "plots/01_oscillator/case2_solution.png"
)

p2_phase = plot(
    x2,
    y2,
    xlabel = "x",
    ylabel = "x'",
    title = "Случай 2 - фазовый портрет",
    label = false,
    lw = 2,
    grid = true
)

savefig(
    p2_phase,
    "plots/01_oscillator/case2_phase.png"
)

p2_combined = plot(
    p2_solution,
    p2_phase;
    layout = (1, 2),
    size = (1200, 450)
)

savefig(
    p2_combined,
    "plots/01_oscillator/case2_combined.png"
)

# ==================================================
# СЛУЧАЙ 3
# С затуханием и внешней силой
#
# x'' + 2x' + x = sin(2t)
#
# Система:
# x' = y
# y' = sin(2t) - 2y - x
# ==================================================

function oscillator3!(du, u, p, t)
    x = u[1]
    y = u[2]

    du[1] = y
    du[2] = sin(2 * t) - 2.0 * y - x
end

prob3 = ODEProblem(
    oscillator3!,
    u0,
    tspan
)

sol3 = solve(
    prob3,
    Tsit5();
    saveat = save_times
)

x3 = [u[1] for u in sol3.u]
y3 = [u[2] for u in sol3.u]

println()
println("="^60)
println("СЛУЧАЙ 3")
println("С затуханием и внешней силой")
println("x'' + 2x' + x = sin(2t)")
println("="^60)

println("x в конце интервала: ", round(x3[end]; digits=6))
println("x' в конце интервала: ", round(y3[end]; digits=6))

p3_solution = plot(
    sol3.t,
    x3,
    xlabel = "t",
    ylabel = "x(t)",
    title = "Случай 3 - решение x(t)",
    label = "x(t)",
    lw = 2,
    grid = true
)

savefig(
    p3_solution,
    "plots/01_oscillator/case3_solution.png"
)

p3_phase = plot(
    x3,
    y3,
    xlabel = "x",
    ylabel = "x'",
    title = "Случай 3 - фазовый портрет",
    label = false,
    lw = 2,
    grid = true
)

savefig(
    p3_phase,
    "plots/01_oscillator/case3_phase.png"
)

p3_combined = plot(
    p3_solution,
    p3_phase;
    layout = (1, 2),
    size = (1200, 450)
)

savefig(
    p3_combined,
    "plots/01_oscillator/case3_combined.png"
)

# --------------------------------------------------
# Завершение
# --------------------------------------------------

println()
println("="^60)
println("РАСЧЁТ ЗАВЕРШЁН")
println("="^60)

println()
println("Сохранены графики:")

println("Случай 1:")
println("  case1_solution.png")
println("  case1_phase.png")
println("  case1_combined.png")

println("Случай 2:")
println("  case2_solution.png")
println("  case2_phase.png")
println("  case2_combined.png")

println("Случай 3:")
println("  case3_solution.png")
println("  case3_phase.png")
println("  case3_combined.png")
