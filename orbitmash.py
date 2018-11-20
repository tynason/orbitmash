#!/usr/bin/python3
# or /usr/bin/env python3
#
# orbitmash.py
# ted nason 2018
# https://github.com/tynason/orbitboy
#___________________________________________________________________#
import matplotlib as mpl
from pylab import *
import scipy.fftpack
import numpy as np
import random
import datetime
import time
import os
import mysql.connector
#___________________________________________________________________#

class Mandel(object):
	"Mandelbrot orbit plotting"
	def __init__(self,doani,doloop,dosave,dodata,chunk,chunksleep,maxiters,minbored,boreme,maxrad,trimend,numbins,numgrads,figclose,figsleep,finalsleep,linesleep,wid,ht,xpos,ypos):
		self.doani=doani;self.doloop=doloop;self.dosave=dosave;self.dodata=dodata;
		self.chunk=chunk;self.chunksleep=chunksleep;self.maxiters=maxiters;
		self.minbored=minbored;self.boreme=boreme;self.maxrad=maxrad;self.trimend=trimend;self.numgrads=numgrads;
		self.figclose=figclose;self.figsleep=figsleep;self.finalsleep=finalsleep;self.linesleep=linesleep;
		self.wid=wid;self.ht=ht;self.xpos=xpos;self.ypos=ypos;

	# https://matplotlib.org/users/customizing.html
	mpl.rcParams['font.size']=8
	mpl.rcParams['axes.labelsize']=8
	mpl.rcParams['xtick.color']='#ffffff'
	mpl.rcParams['ytick.color']='#ffffff'
	mpl.rcParams['lines.linewidth']=0.4

	def getcolor(self):
		brite=0
		while True:
			r=np.random.uniform(0.4,0.8)
			g=np.random.uniform(0.0,0.2)
			b=np.random.uniform(0.4,0.8)
			# not too dark
			if r>0.6 or b>0.6: break
		print('rgbfore: ',int(r*256),int(g*256),int(b*256))
		return (r,g,b)

	def lumengen(self,fore,back,grads): # get colors between start (back) and end (fore) colors
		rangeRGB=[x1-x2 for (x1,x2) in zip(fore,back)] # the range of RGB[0-1] to be covered
		segRGB=list(map(lambda x: x/grads,rangeRGB)) # the amount to decrement each element RGB
		R=np.zeros((grads,));G=np.zeros((grads,));B=np.zeros((grads,)) # start w/fore and decrement to back
		for nn in range(self.numgrads): R[nn]=fore[0]-nn*segRGB[0];G[nn]=fore[1]-nn*segRGB[1];B[nn]=fore[2]-nn*segRGB[2]
		return list(zip(R,G,B))

	def getpq(self):
		print('searching...')
		boredcount=0
		while True:
			boredcount+=1
			params=[]
			flip=x=np.random.random_integers(1,5)
			if flip==1: # 5 period repulsive
				p=np.random.uniform(-0.45,-0.6);q=np.random.uniform(0.5,0.54)
			elif flip==2: # 5 period attractive
				p=np.random.uniform(0.355,0.360);q=np.random.uniform(0.315,0.400)

			elif flip==3: # elephant valley
				p=np.random.uniform(0.2,0.4);q=np.random.uniform(0,0.4)
			elif flip==4: # seahorse valley
				p=np.random.uniform(-0.8,-0.7);q=np.random.uniform(0,0.1)
			elif flip==5: # MAIN CARDIOID
				theta=np.random.uniform(0,2*math.pi);rrr=np.random.uniform(0,0.05)
				p=(0.5+rrr)*cos(theta)-cos(2*theta)/4;q=(0.5+rrr)*sin(theta)-sin(2*theta)/4
			elif flip==6: # BULB AT (-1.0)
				theta=np.random.uniform(0,2*math.pi);rrr=np.random.uniform(0,0.05)
				p=(0.25+rrr)*cos(theta)-1;q=(0.25+rrr)*sin(theta)
			#___________________________________________________________________#


			k=0
			c=complex(p,q)
			z=complex(0.0,0.0)

			while abs(z)<self.maxrad and k<maxiters:
				k+=1;z=z*z+c

			if self.boreme:
				if k>=self.minbored: # it's boring but over the minbored, so stop looking
					print('BORING #',boredcount,'kappa=',k,'\tp=',p,'\tq=',q)
					break
				else: continue # it's less then minbored so keep looking
			else:
				if k>self.minbored and k<maxiters: # we found a non-boring one
					break
				else: continue # keep looking

		# we return a p,q that is not boring
		return [p,q,k]

	def datagen(self,p,q):
		zdata=[]
		k=0
		c=complex(p,q)
		z=complex(0.0,0.0)
		zdata.append(z)

		while abs(z)<self.maxrad and k<maxiters:
			k+=1;z=z*z+c
			zdata.append(z)

		for i in range(trimend):
			zdata.pop()

		#print(zdata)
		return zdata

	def refreshall(self):
		for ax in fig.axes:
			currcolor = ax.patch.get_facecolor();currtitle = ax.get_title(loc='left')
			currxlabel=ax.get_xlabel();currylabel=ax.get_ylabel()
			currxlim=ax.get_xlim();currylim=ax.get_ylim()
			ax.clear();ax.patch.set_facecolor(currcolor)
			ax.set_title(currtitle,loc='left',color=mybritegrn)
			ax.grid(True);ax.patch.set_alpha(1.0)
			ax.set_xlabel(currxlabel,color=mybritegrn)
			ax.set_ylabel(currylabel,color=mybritegrn)
			#ax.set_xlim(currxlim);ax.set_ylim(currylim)

	def refreshone(self,n):
		ax=fig.axes[n]
		currcolor = ax.patch.get_facecolor();currtitle = ax.get_title(loc='left')
		currxlabel=ax.get_xlabel();currylabel=ax.get_ylabel()
		currxlim=ax.get_xlim();currylim=ax.get_ylim()
		ax.clear();ax.patch.set_facecolor(currcolor)
		ax.set_title(currtitle,loc='left',color=mybritegrn)
		ax.grid(True);ax.patch.set_alpha(1.0)
		ax.set_xlabel(currxlabel,color=mybritegrn)
		ax.set_ylabel(currylabel,color=mybritegrn)
		ax.set_xlim(currxlim);ax.set_ylim(currylim)

#___________________________________________________________________#
start_time=time.time()

mygunmet='#113344';myblue='#11aacc';mydkblue='#0000cc'
myturq='#00ffff';myturq2='#11bbbb';myteal='#00ffcc';myteal2='#00ccaa'
mygreen='#44bb44';mybritegrn='#00ff66';myfuscia='#ff00ff';myfuscia2='#dd0099'
mypurp='#ff00cc';mypurp2='#9933ff';myred='#cc0000';myorange='#ffaa00'
myorange2='#ff6600';myyell='#ffff00';myyell2='#ffcc00'

mygunmet2='#052529'
rgbback=mygunmet2
#___________________________________________________________________#

doani=False 		# animate the orbit, rad, and ang plots
doloop=True		# loop thru random orbits, not just one
dosave=False 	# save params to DB or file, and save png

dodata=False
chunk=40
chunksleep=0

maxiters=9000	# max iterations
minbored=2000	# minimum non-boring orbit iterations

boreme=False 	#  False to pick a long orbit which escapes at the end
maxrad=2.0		# defines the escape criterion

trimend=4		# omits the final few iterations from some of the plots
numbins=120		# no. of bins in the histograms
numgrads=12		# how many slices to plot the animated orbit

figclose=True 	# close after plotting
figsleep=4		# various sleep intervals
finalsleep=3
linesleep=0

wid=1800;ht=1100;xpos=30;ypos=30
#___________________________________________________________________#

mymand=Mandel(doani,doloop,dosave,dodata,chunk,chunksleep,maxiters,minbored,boreme,maxrad,trimend,numbins,numgrads,figclose,figsleep,finalsleep,linesleep,wid,ht,xpos,ypos)
print(mymand.__doc__)

#___________________________________________________________________#

fig = plt.figure()

# initial subplot settings
ax0 = plt.subplot2grid((3,4),(0,0))
ax1 = plt.subplot2grid((3,4),(0,1))
ax2 = plt.subplot2grid((3,4),(0,2))
ax3 = plt.subplot2grid((3,4),(1,0))
ax4 = plt.subplot2grid((3,4),(1,1))
ax5 = plt.subplot2grid((3,4),(1,2))
ax6 = plt.subplot2grid((3,4),(2,0))
ax7 = plt.subplot2grid((3,4),(2,1))
ax8 = plt.subplot2grid((3,4),(2,2))

ax9 = plt.subplot2grid((3,4),(0,3))
ax10 = plt.subplot2grid((3,4),(1,3))
ax11 = plt.subplot2grid((3,4),(2,3))

for ax in fig.axes:
	ax.grid(True)
	ax.patch.set_facecolor(rgbback)
	ax.patch.set_alpha(1.0)

win = plt.gcf().canvas.manager.window
fig.canvas.manager.window.wm_geometry('%dx%d%+d%+d' % (wid,ht,xpos,ypos))
fig.patch.set_facecolor(rgbback)
plt.tight_layout(pad=0.1,w_pad=0.1,h_pad=0.5)
fig.subplots_adjust(top=0.935)
#___________________________________________________________________#


"""
mymand.getcolor
mymand.lumengen
"""


while True:

	result=mymand.getpq()
	p=result[0]
	q=result[1]
	kappa=result[2]
	c=complex(p,q)
	print('\np,q:\t',c.real,c.imag);print('\nc:\t',c) #print('\nzdata:\t',zdata)

	delta=0.005  # generate 8 more values of p,q
	p0=p-delta;q0=q+delta;c0=complex(p0,q0)
	p1=p;q1=q+delta;c1=complex(p1,q1)
	p2=p+delta;q2=q+delta;c2=complex(p2,q2)
	p3=p-delta;q3=q;c3=complex(p3,q3)
	c4=c
	p5=p+delta;q5=q;c5=complex(p5,q5)
	p6=p-delta;q6=q-delta;c6=complex(p6,q6)
	p7=p;q7=q-delta;c7=complex(p7,q7)
	p8=p+delta;q8=q-delta;c8=complex(p8,q8)

	carray=np.array([c0,c1,c2,c3,c4,c5,c6,c7,c8],dtype=complex)

	print('\nSELECTED non-boring p,q:')
	print('\tp='+str(p));print('\tq='+str(q));print('\tboreme:',boreme);print('\tminbored:',minbored)
	print('\tmaxiters',maxiters);print('\tnumgrads',numgrads);print('\ttrimend:',trimend)

	mandfunc='z*z+c'
	mytitle='Mandelbrot ' +mandfunc+ '  p= '+str(p)+'  q= '+str(q)+'  kappa= '+str(kappa)+'  maxrad= '+str(maxrad)
	fig.suptitle(mytitle,fontsize=9,color='#00ff00')

	mymand.refreshall()
	ax0.set_title('x-y orbit plot  trimend='+str(trimend),loc='left',color=mybritegrn)

	for i in range(0,9):
		#print('\nplotting axes no.: \t',i)

		zdata=mymand.datagen(carray[i].real,carray[i].imag)

		axes=i
		ax=fig.axes[axes]
		X = [x.real for x in zdata]
		Y = [x.imag for x in zdata]
		for i in range(trimend):
			X.pop()
			Y.pop()

		ax.plot(X,Y,color=mybritegrn)

		ax9.clear();ax9.patch.set_facecolor(mygunmet2)
		ax9.hist(X,bins=numbins,normed=True,color=myyell)

		ax10.clear();ax10.patch.set_facecolor(mygunmet2)
		ax10.plot(X,color=mybritegrn)

		ax11.clear();ax11.patch.set_facecolor(mygunmet2)


		# do some fourier stuff
		data=X
		kappa=len(X)
		datalbl='rad'
		maxfreq=int(kappa/15)

		Fs = 300  # sampling rate - HIGHER = SHARPER PEAK
		Ts = 1.0/Fs # sampling interval
		n=kappa-2 # length of the signal
		kk = np.arange(kappa-2)
		T = kappa/Fs
		frq2 = kk/T # two sides frequency range
		frq1 = frq2[1:int(n/2)] # one side frequency slice
		Y2 = 2*np.fft.fft(data)/n # fft computing and norm
		Y1 = Y2[1:int(n/2)]

		ax11.set_title('FT('+datalbl+') vs period   Fsamp=300',loc='left',color=mybritegrn)
		ax11.set_xlabel('period',color=mybritegrn)
		ax11.set_ylabel('FT('+datalbl+')',color=mybritegrn)
		ax11.plot(frq1[1:maxfreq],abs(Y1[1:maxfreq]),color=myorange)

		plt.show(block=False);fig.canvas.draw()
	#___________________________________________________________________#

	time.time();elapsed=time.time()-start_time
	timerpt='\telapsed: '+str(elapsed);print(timerpt)
	time.sleep(figsleep)

	if not doloop: break

if figclose:
	time.sleep(finalsleep)
	plt.close()
else:
	plt.show(block=True)
































#___________________________________________________________________#
