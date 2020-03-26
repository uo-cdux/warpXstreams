import sys
import os
import os.path
import numpy
import pandas
import itertools
import toyplot
import toyplot.pdf

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

canvas = toyplot.Canvas('4in','2.5in')
axes = canvas.cartesian(xlabel = 'A.I. (FLOP/Byte)',
                        ylabel='Performance (GFLOP/sec)', margin=(50, 50))
x = data['intensity']
y = data['flop']
mark = axes.scatterplot(x, y)
toyplot.pdf.render(canvas, "%s.pdf" %filename)
