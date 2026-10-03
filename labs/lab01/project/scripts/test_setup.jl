using DrWatson

@quickactivate "project"

println("Проект активирован: ", projectdir())

packages = [
    "DrWatson",
    "DifferentialEquations",
    "Plots",
    "DataFrames",
    "CSV",
    "JLD2",
    "Literate",
    "IJulia",
    "BenchmarkTools",
    "Quarto"
]

println()
println("Проверка пакетов:")

for pkg in packages
    try
        eval(Meta.parse("using $pkg"))
        println(" ✓ $pkg")
    catch e
        println(" ✗ $pkg: Ошибка загрузки")
    end
end

println()
println("Структура проекта:")
println(" Корень: ", projectdir())
println(" Данные: ", datadir())
println(" Скрипты: ", scriptsdir())
println(" Графики: ", plotsdir())
