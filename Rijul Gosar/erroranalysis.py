# -*- coding: utf-8 -*-
"""
Created on Mon Feb 01 18:31:05 2016

@author: Eeshani Godbole
"""

import scipy,numpy
import scipy.optimize, scipy.stats
import numpy.random
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels
import statsmodels.stats
import statsmodels.stats.stattools as stools

plt.style.use("ggplot")


def formataxis(ax):
    ax.xaxis.label.set_fontname('Georgia')
    ax.xaxis.label.set_fontsize(12)
    ax.yaxis.label.set_fontname('Georgia')
    ax.yaxis.label.set_fontsize(12)
    ax.title.set_fontname('Georgia')
    ax.title.set_fontsize(12)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(8)
        
def get_stderr_fit(f,Xdata,popt,pcov,dict_data):
    Y=f(Xdata,popt,dict_data)
    listdY=[]
    for i in xrange(len(popt)):
        p=popt[i]
        dp=abs(p)/1e6 + 1e-20
        popt[i] += dp
        Yi=f(Xdata,popt,dict_data)
        dY=(Yi-Y)/dp
        listdY.append(dY)
        popt[i] -=dp
    listdY=scipy.array(listdY)
    left=scipy.dot(listdY.T,pcov)
    right=scipy.dot(left,listdY)
    sigma2y=right.diagonal()
    mean_sigma2y=scipy.mean(right.diagonal())
    M=Xdata.shape[1]
    N=len(popt)
    avg_stddev_data=scipy.sqrt(M*mean_sigma2y/N)
    sigmay=scipy.sqrt(sigma2y)
    return sigmay,avg_stddev_data
    popt=numpy.array([0.0]*11)
    
def fitdata(f,Xdata,Ydata,Errdata,pguess,dict_data,ax=False,ax2=False):
    
    def error(p,Xdata,Ydata,Errdata,dict_data):
        Y=f(Xdata,p,dicct_data)
        residuals=(Y-Ydata)/Errdata
        return residuals
    res=scipy.optimize.leastsq(error,pguess.args=(Xdata,Ydata,Errdata,dict_data),full_output=1)
    (popt,pcov,infodict,errmsg,ier)=res
    perr=scipy.sqrt(scipy.diag(pcov))
    
    M=len(Ydata)
    N=len(popt)
    
    Y=f(Xdata,popt,dict_data)
    residuals=(Y-Ydata)/Errdata
    meanY=scipy.mean(Ydata)
    squares=(Y-meanY)/Errdata
    squaresT=(Ydata-meanY)/Errdata
    
    SSM=sum(squares**2)
    SSE=sum(residuals**2)
    SST=sum(squaresT**2)
    
    DFM=N-1
    DFE=M-N
    DFT=M-1
    
    MSM=SSM/DFM
    MSE=SSE/DFE
    MST=SST/DFT
    
    R2=SSM/SST
    R2_adj=1-(1-R2)*(M-1)/(M-N-1)
    
    t_stat=popt/perr
    t_stat=t_stat.real
    p_p=1.0-scipy.stats.t.cdf(t_stat,DFE)
    z=scipy.stats.t(M-N).ppf(0.95)
    p95=perr*z
    
    chisquared=sum(residuals**2)
    degfreeedom=M-N
    chisqured_red=chisquared/degfreedom
    p_chi2=1.0-scipy.stats.chi2.cdf(chisquared,degfreedom)
    chisquare=(p_chi2,chisquared,chisquared_red,degfreedom,R2,R2_adj)
    
    w,p_shapiro=scipy.stats.shapiro(residuals)
    mean_res=scipy.mean(residuals)
    stddev_res=scipy.sqrt(scipy.var(residuals))
    t_res=mean_res/stddev_res
    p_res=1.0-scipy.stats.t.cdf(t_res,M-1)
    
    F=MSM/MSE
    p_F=1.0-scipy.stats.f.cdf(F,DFM,DFE)
    
    dw=stools.durbin_watson(residuals)
    resanal=(p_shapiro,w,mean_res,p_res,F,p_F,dw)
    
    if ax:
        formataxis(ax)
        ax.plotdata(Ydata,Y,'ro')
        ax.errorbar(Ydata,Y,yerr=Errdata,fmt='.')
        Ymin,Ymax=min((min(Y),min(Ydata))),max((max(Y),max(Ydata)))
        ax.plot([Ymin,Ymax],[Ymin,Ymax],'b')
        
        ax.xaxis.label.set_text('Data')
        ax.yaxis.label.set_text('Fitted')
        
        sigmay,avg_stddev_data=get_stderr_fit(f,Xdata,popt,pcov,dict_data)
        Yplus=Y+sigmay
        Yminus=Y-sigmay
        ax.plot(Y,Yplus,'c',alpha=0.6,linestyle='--',linewidth=0.5)
        ax.plot(Y,Yminus,'c',alpha=0.6,linestyle='--',linewidth=0.5)
        ax.fill_between(Y,Yminus,Yplus,facecolor='cyan',alpha=0.5)
        titletext ='parity plot for fit.\n'
        titletext +=r'$r^2$=%5.2f,$r^2_{adj}$=%5.2f,'
        titletext +='$\sigma_{exp}$=%5.2f,$\chi^2_{\nu}$=%5.2f,$p_{\chi^2}$=%5.2f,'
        titletext +='$\sigma_{err}^{reg}$=%5.2f'
        
        ax.title.set_text(titletext%(R2,R2_adj,avg_stddev_data,chisquared_red,p_chi2,stderr_reg))
        ax.figure.canvas.draw()
        
    if ax2:
        formataxis(ax2)
        ax2.plot(Y,residuals,'ro')
        
        ax2.xaxis.label.set_text('Fitted Data')
        ax2.yaxis.label.set_text('Residuals')
        
        titletext ='Analysis of Residuals\n'
        titletext +=r'mean=%5.2f,$p_{res}$=%5.2f,$p_{shapiro}$=%5.2f, $Durbin-Watson=%2.1f'
        titletext +='\n F=%5.2f, $p_F$=%3.2e'
        ax2.title.set_text(titletext%(mean_res,p_res,p_shapiro,dw,F,p_F))
        
        ax2.figure.canvas.draw()
        
    return popt,pcov,perr,p95,p_p,chisquare,resanal

def import_data(xlfile,sheetname):
    df=pd.read_excel(xlfile,sheetname=sheetname)
    return df 

def prepare_data(df, Criterion, Predictors, Error=False):
    Y=scipy.array(df[Criterion])
    if Error:
        Errdata=scipy.array(df[Error])
    else:
        Errdata=scipy.ones(len(Y))
    Xdata=[]
    for X in Predictors:
        X=list(df[X])
        Xdata.append(X)
    Xdata=scipy.array(Xdata)
    return Xdata, Y, Errdata
    
if  _name_=="_main_":
    fig=plt.figure()
    ax=fig.add_subplot(111)
    fig.show()
    
    fig2=plt.figure()
    ax2=fig2.add_subplot(111)
    fig2.show()
    
    
    def f(X,p,dict_data):
        a=dict_data['a']
        b=dict_data['b']
        (x,y,z)=X
        Y=p[0]+p[1]*X**2+p[2]*y+p[3]*z
        return Y
        
    df=import_data('SynthData.xlsx','Data')
    Xdata, Ydata,Errdata=prepare_data(df, 'Ydata',('x','y','z'),Error='err')
    
    N=4
    pguess=N*[0.0]
    
    popt, pcov, perr, p95, p_p, chisquare, resanal=fitdata(f,Xdata,Ydata,Errdata,pguess,dict_data,ax=ax,ax2=ax2)
    
    