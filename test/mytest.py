from dspyields.dspyields import maximum_yield

kegg_pathway = ['R01788', 'R02737', 'R02780', 'R02628', 'R00724']
organism = "ecoli"
model = "iJO1366.xml"

maximum_yield(kegg_pathway, organism, model)