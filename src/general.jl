import Dates

timestamp() = Dates.format(Dates.now(), "yymmdd_HMS")

function savefig(basename, fig)
    save(plotsdir(basename * ".png"), fig; px_per_unit=2.0)
end
