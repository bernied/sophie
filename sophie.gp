
h / help				flag	"  display this help and exit"
V / version				flag	"  output version information and exit"
p / print-path			flag    "  Print tour for lowest time"
"Options controlling simulated annealing heuristics:"
b / brute-force			flag	"  Use brute force to solve"
i /initial-temp			float	"  Set starting temp for simulated annealing"
t / steps-per-temp		int		"  Steps per temperature change"
s / cooling-steps		int		"  Numer of cooling steps"
f / cooling-fraction	float	"  Cooling fraction for exponent"
k / K					float	"  K factor"
r / runs				int		"  Number of runs"
e / feed-back           flag    "  Feed successive runs back into the next run"
v / max-verts			int		"  How many verticies before using heuristic"
z / randomize			flag	"  Randomize starting heuristic path"

#usage_begin
Usage: __PROGRAM_NAME__ [OPTION]... [FILE]

__GLOSSARY_GNU__

#usage_end

