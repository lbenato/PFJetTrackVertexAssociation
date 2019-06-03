variable = {}

var_template = {
    "EventNumber": {
      "title" : "event number",
      "nbins" : 10000000,
      "min" : 0,
      "max" : 1.e7,
      "log" : False,
    },
    "LumiNumber": {
      "title" : "lumisection number",
      "nbins" : 2000,
      "min" : 0,
      "max" : 2000,
      "log" : False,
    },
    "RunNumber": {
      "title" : "run number",
      "nbins" : 7000,
      "min" : 254000,
      "max" : 261000,
      "log" : False,
    },
    "nPV": {
      "title" : "number of reconstructed Primary Vertices",
      "nbins" : 50,
      "min" : -0.5,
      "max" : 49.5,
      "log" : False,
    },
    "nJets": {
      "title" : "number of jets",
      "nbins" : 17,
      "min" : -0.5,
      "max" : 16.5,
      "log" : True,
    },
}



for n, v in var_template.iteritems():
    if '[N]' in n:
        for i in range(0, 5):
            ni = n.replace('[N]', "%d" % i)
            variable[ni] = v.copy()
            variable[ni]['title'] = variable[ni]['title'].replace('[N]', "%d" % i)
    else:
        variable[n] = v
