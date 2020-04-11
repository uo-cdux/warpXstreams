import sys
import os
import os.path
import numpy
import pandas
import itertools
import toyplot
import toyplot.pdf

colors = ['#00429d', '#204fa3', '#325ca9', '#4169ae', '#4e77b2',
          '#5b84b7', '#6892ba', '#74a0bd', '#82aebf', '#90bcbf',
          '#9fcbbf', '#afd9bb', '#c1e7b5', '#d7f4a7', '#fbff7c']
colors.reverse()

def GetTimingsFromFile(data):
  toWrite = []
  FLOPS = []
  INTENSITY = []
  PARAMS = ""
  for line in data:
    if "PARAMS" in line :
      if len(PARAMS) is 0:
        PARAMS = (line.split(":")[1]).strip()
        FLOPS = []
        INTENSITY = []
      else :
        flops = numpy.mean(FLOPS)
        intensity = numpy.mean(INTENSITY)
        toWrite.append("%s,%f,%f\n" %(PARAMS,flops,intensity))
        PARAMS = (line.split(":")[1]).strip()
        FLOPS = []
        INTENSITY = []
    else :
      if len(line) is not 1 and len(line) is not 0:
        linedata = line.split(",") #Separate GFLOPs and FLOP/Dype
        flops = float((linedata[0].split(":")[1]).strip())
        intensity = float((linedata[1].split(":")[1]).strip())
        FLOPS.append(flops)
        INTENSITY.append(intensity)
  flops = numpy.mean(FLOPS)
  intensity = numpy.mean(INTENSITY)
  toWrite.append("%s,%f,%f\n" %(PARAMS,flops,intensity))
  return toWrite

outfile = "data.csv"
output = open(outfile, "w")
header = "params,flop,intensity\n"
output.write(header);

if len(sys.argv) < 2 :
  print "Input file not provided"
  sys.exit(-1)


filename = sys.argv[1]
if os.path.isfile(filename) and os.access(filename, os.R_OK):
  data = open(filename, "r")
  toWrite = GetTimingsFromFile(data)
  for line in toWrite :
    output.write(line)
  output.close()
  data.close()

data = pandas.read_csv(outfile)
#forplots = data.pivot_table(values='time', columns=['datasets','algorithm'], index='parenv')
print data

markers = []
indices = range(0, 15)
canvas = toyplot.Canvas('4in','2.5in')
axes = canvas.cartesian(xlabel = 'A.I. (FLOP/Byte)',
                        ylabel='Performance (GFLOP/sec)', margin=(20, 100, 50, 50))
for (x, y, index) in zip(data['intensity'], data['flop'], indices) :
  mark = axes.scatterplot(x, y, color=colors[index], size=8, opacity=1.0)
  markers.append(mark)
legend = zip(data['params'], markers)
canvas.legend(legend, corner=("right", 55, 50, 150))
toyplot.pdf.render(canvas, "%s.pdf" %filename)
