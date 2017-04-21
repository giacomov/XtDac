#Author: Giacomo Vianello (giacomov@stanford.edu)

import urllib, base64
import StringIO
import matplotlib.pyplot as plt
import time
import numpy

htmlCode1 = '''
<html>
<head>
<style>
.CSSTableGenerator {
	margin:0px;padding:0px;
	width:100%;
	box-shadow: 10px 10px 5px #888888;
	border:1px solid #000000;
	
	-moz-border-radius-bottomleft:10px;
	-webkit-border-bottom-left-radius:10px;
	border-bottom-left-radius:10px;
	
	-moz-border-radius-bottomright:10px;
	-webkit-border-bottom-right-radius:10px;
	border-bottom-right-radius:10px;
	
	-moz-border-radius-topright:10px;
	-webkit-border-top-right-radius:10px;
	border-top-right-radius:10px;
	
	-moz-border-radius-topleft:10px;
	-webkit-border-top-left-radius:10px;
	border-top-left-radius:10px;
}.CSSTableGenerator table{
    border-collapse: collapse;
        border-spacing: 0;
	width:100%;
	height:100%;
	margin:0px;padding:0px;
}.CSSTableGenerator tr:last-child td:last-child {
	-moz-border-radius-bottomright:10px;
	-webkit-border-bottom-right-radius:10px;
	border-bottom-right-radius:10px;
}
.CSSTableGenerator table tr:first-child td:first-child {
	-moz-border-radius-topleft:10px;
	-webkit-border-top-left-radius:10px;
	border-top-left-radius:10px;
}
.CSSTableGenerator table tr:first-child td:last-child {
	-moz-border-radius-topright:10px;
	-webkit-border-top-right-radius:10px;
	border-top-right-radius:10px;
}.CSSTableGenerator tr:last-child td:first-child{
	-moz-border-radius-bottomleft:10px;
	-webkit-border-bottom-left-radius:10px;
	border-bottom-left-radius:10px;
}.CSSTableGenerator tr:hover td{
	
}
.CSSTableGenerator tr:nth-child(odd){ background-color:#fcc694; }
.CSSTableGenerator tr:nth-child(even)    { background-color:#ffffff; }.CSSTableGenerator td{
	vertical-align:middle;
	
	
	border:1px solid #000000;
	border-width:0px 1px 1px 0px;
	text-align:center;
	padding:7px;
	font-size:12px;
	font-family:Arial;
	font-weight:normal;
	color:#000000;
}.CSSTableGenerator tr:last-child td{
	border-width:0px 1px 0px 0px;
}.CSSTableGenerator tr td:last-child{
	border-width:0px 0px 1px 0px;
}.CSSTableGenerator tr:last-child td:last-child{
	border-width:0px 0px 0px 0px;
}
.CSSTableGenerator tr:first-child td{
		background:-o-linear-gradient(bottom, #ff7f00 5%, #bf5f00 100%);	background:-webkit-gradient( linear, left top, left bottom, color-stop(0.05, #ff7f00), color-stop(1, #bf5f00) );
	background:-moz-linear-gradient( center top, #ff7f00 5%, #bf5f00 100% );
	filter:progid:DXImageTransform.Microsoft.gradient(startColorstr="#ff7f00", endColorstr="#bf5f00");	background: -o-linear-gradient(top,#ff7f00,bf5f00);

	background-color:#ff7f00;
	border:0px solid #000000;
	text-align:center;
	border-width:0px 0px 1px 1px;
	font-size:14px;
	font-family:Arial;
	font-weight:bold;
	color:#ffffff;
}
.CSSTableGenerator tr:first-child:hover td{
	background:-o-linear-gradient(bottom, #ff7f00 5%, #bf5f00 100%);	background:-webkit-gradient( linear, left top, left bottom, color-stop(0.05, #ff7f00), color-stop(1, #bf5f00) );
	background:-moz-linear-gradient( center top, #ff7f00 5%, #bf5f00 100% );
	filter:progid:DXImageTransform.Microsoft.gradient(startColorstr="#ff7f00", endColorstr="#bf5f00");	background: -o-linear-gradient(top,#ff7f00,bf5f00);

	background-color:#ff7f00;
}
.CSSTableGenerator tr:first-child td:first-child{
	border-width:0px 0px 1px 0px;
}
.CSSTableGenerator tr:first-child td:last-child{
	border-width:0px 0px 1px 1px;
}
</style>
</head>
<body>
#JOB_INFO#
#OBS_INFO#
#HWU_LC#
<h2>Excesses</h2>
<div class="CSSTableGenerator" >
<table>
<tr>
<td>
Region info
</td>
<td>
Block edges
</td>
<td>
Finding map
</td>
<td>
Stamps
</td>
<td>
Light curve
</td>
</tr>
'''
dataHtml = '''
<tr><td>#INFO#</td><td>#INTERVALS#</td><td><img src=#MAP#/></td><td><img src=#STAMPS#/></td><td><img src=#LC#/></td></tr>
'''

htmlCode2 = '''
</table>
</div>

</body>
</html>
'''

class Summary(object):
  def __init__(self):
    self.htmlCode1            = str(htmlCode1)
    self.dataHtml             = ''
    self.htmlCode2            = str(htmlCode2)
  pass
  
  def addJobInfo(self,runtime,ncpu,typeIprob):
    jobInfo                   = '<h1>EXTraS Divide And Conquer: summary of results</h1>'
    jobInfo                  += '<h2>Job info</h2>'
    jobInfo                  += 'Date: %s<br/>' % time.asctime()
    jobInfo                  += 'Runtime: %.2f<br/>' %runtime
    jobInfo                  += 'Number of cores: %s<br/>' % ncpu
    jobInfo                  += 'Type I error prob.: %s<br/>' % typeIprob
    self.htmlCode1            = self.htmlCode1.replace("#JOB_INFO#",jobInfo)
  
  def addObsInfo(self,filename,instrumentName,tstart,tstop,nEvents):
    obsInfo                   = "<h2>Observation info</h2>"
    obsInfo                  += "File name: %s<br/>" % filename
    obsInfo                  += "Instrument: %s<br/>" % instrumentName
    obsInfo                  += "Start time: %s<br/>" % tstart
    obsInfo                  += "Stop time: %s<br/>" % tstop
    obsInfo                  += "Duration: %.3f<br/>" % (tstop-tstart)
    obsInfo                  += "Number of events: %s<br/>" % (nEvents)
    self.htmlCode1            = self.htmlCode1.replace("#OBS_INFO#",obsInfo)
  
  def addWholeHWUlightCurve(self,time,tstart,tstop,**kwargs):
    NBins                     = int(numpy.ceil((tstop-tstart)/50.0))
        
    (counts,edges)            = numpy.histogram(time,NBins)
    leftEdges                 = edges[:-1]
    width                     = ((tstop-tstart)/NBins)
    rates                     = counts / width
    
    fig                       = plt.figure(**kwargs)
    sub                       = fig.add_subplot(111)

    rr = numpy.append(rates, 0)

    sub.step(edges - edges[0], rr, where='post')
    sub.set_xlabel("Mission Elapsed Time (s)")
    sub.set_ylabel("Rate")
    sub.set_title("Light curve for the whole HWU")
    
    hwulc                     = StringIO.StringIO()
    fig.savefig(hwulc,tight_layout=True,format='png')
    plt.close(fig)
    hwulc.seek(0)
    
    value                     = urllib.quote(base64.b64encode(hwulc.buf))
    
    hwu_lc                    = "<h2>Light curve for the whole Hardware Unit</h2>"
    hwu_lc                   += "<img src='data:image/jpg;base64," + value + "'>"
    self.htmlCode1            = self.htmlCode1.replace("#HWU_LC#",hwu_lc)
    
  def addResults(self,interestingRegions):
    #Obtain the pretty pictures from the InterestingRegion instances
    #and add them to the HTML code
    for i,interestingRegion in enumerate(interestingRegions):
      fig                      = interestingRegion.getLightCurve()
      lc                       = StringIO.StringIO()
      fig.savefig(lc,tight_layout=True,format='png')
      plt.close(fig)
      lc.seek(0)
      
      fig                      = interestingRegion.getStamps()
      stamps                   = StringIO.StringIO()
      fig.savefig(stamps,tight_layout=True,format='png')
      plt.close(fig)
      stamps.seek(0)
      
      fig                      = interestingRegion.getFindingMap()
      cmap                     = StringIO.StringIO()
      fig.savefig(cmap,tight_layout=True,format='png')
      plt.close(fig)
      cmap.seek(0)
      
      ra,dec                   = interestingRegion.getCenterSkyCoordinates()
      x,y                      = (interestingRegion.box.c1, interestingRegion.box.c2)
      w,h                      = interestingRegion.getWidthAndHeight()
      intervalsEdges           = interestingRegion.getIntervalsEdges()
      interestingRegion.clearMemory()
      
      self.addOneResult(cmap,stamps,lc,ra,dec,x,y,w,h,intervalsEdges)
    pass
  pass
  
  def addOneResult(self,cmap,stamps,lc,ra,dec,x,y,width,height,intervalsEdges):
    cmapdata                  = urllib.quote(base64.b64encode(cmap.buf))
    stampsdata                = urllib.quote(base64.b64encode(stamps.buf))
    lcdata                    = urllib.quote(base64.b64encode(lc.buf))
    
    thisDataHtml              = str(dataHtml)
    for key,value in {'#MAP#':cmapdata,'#STAMPS#': stampsdata, '#LC#': lcdata}.iteritems():
      thisDataHtml            = thisDataHtml.replace(key,"'data:image/jpg;base64," + value + "'")
    pass
    
    #Now build the INFO column
    info                      = "Center (R.A., Dec.): (%.5f,%.5f)<br/>" %(ra,dec)
    info                     += "Width x height     : %.2f x %.2f px<br/>" %(width,height)
    info                     += "Edge (X,Y)       : (%.2f,%.2f)<br/>" %(x,y)
    thisDataHtml              = thisDataHtml.replace("#INFO#",info)
    
    #Now build the INTERVALS column
    intervals                 = ("%s" % 
                                 ("".join(map(lambda x:"%.3f<br/>" % x,intervalsEdges))))
    thisDataHtml              = thisDataHtml.replace("#INTERVALS#",intervals)
    
    self.dataHtml            += thisDataHtml
    
  pass
  
  def write(self,filename):
    with open(filename,"w+") as f:
      f.write(self.htmlCode1)
      f.write(self.dataHtml)
      f.write(self.htmlCode2)
    pass
  pass
