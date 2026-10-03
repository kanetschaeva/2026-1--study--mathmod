# # Экспоненциальный рост
#
# **Цель:** исследовать решение уравнения du/dt = αu.
#
# Экспоненциальный рост - это процесс увеличения величины,
# при котором скорость роста в каждый момент времени
# пропорциональна текущему значению этой величины.

using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots
using DataFrames
using JLD2

script_name = splitext(basename(PROGRAM_FILE))[1]

mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

# ## Определение модели

function exponential_growth!(du, u, p, t)
    α = p
    du[1] = α * u[1]
end

# ## Первый запуск

u0 = [1.0]
α = 0.3
tspan = (0.0, 10.0)

prob = ODEProblem(
    exponential_growth!,
    u0,
    tspan,
    α
)

sol = solve(
    prob,
    Tsit5(),
    saveat=0.1
)

# ## Визуализация

plot(
    sol,
    label="u(t)",
    xlabel="Время t",
    ylabel="Популяция u",
    title="Экспоненциальный рост (α = $α)",
    lw=2,
    legend=:topleft
)

savefig(
    plotsdir(
        script_name,
        "exponential_growth_α=$α.png"
    )
)

# ## Анализ

df = DataFrame(
    t=sol.t,
    u=first.(sol.u)
)

println("Первые 5 строк результатов:")
println(first(df, 5))

u_final = last(sol.u)[1]
doubling_time = log(2) / α

println()
println(
    "Финальное значение u(t): ",
    round(u_final; digits=3)
)

println(
    "Аналитическое время удвоения: ",
    round(doubling_time; digits=2)
)

# ## Сохранение результатов

@save datadir(
    script_name,
    "all_results.jld2"
) df
