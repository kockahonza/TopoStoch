import Dates

timestamp() = Dates.format(Dates.now(), "yymmdd_HMS")

function savefig(subdir, basename, fig)
    save(plotsdir(subdir, basename * ".png"), fig; px_per_unit=2.0)
end
savefig(basename, fig) = savefig("", basename, fig)
