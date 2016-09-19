#Matplotlib configuration for xtdac

rcParams = {}

rcParams['axes.facecolor'] = 'white'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'x-large'
rcParams['figure.dpi']     = 50
rcParams['figure.figsize'] = [5,5]
rcParams['font.size'] = 8
rcParams['interactive'] = False
rcParams['savefig.dpi'] = 50
rcParams['xtick.labelsize'] = 'x-large'
rcParams['ytick.labelsize'] = 'x-large'

def getConfig():
  return rcParams
