

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
    '''
    popt,pcov,perr,chisquare,shapiro=get_stderr_fit(f,xdata,popt,pcov,dict_data)
    
    f=function f(x,p,dict_data)
    Xdata=array like object (k,M) shaped array for data with k predictors
    e.g. if X=(x1,x2,x3) then X=(x1,x2,x3) where x1 is a vector of x1 etc
    ydata=array like object of length M
    popt=vector of length N of the  optimized parameters
    pcov=covariance matrix of the fit
    dict_data=dictionary containing other data necessary for f
    
    returns
    sigmay: array like object of length M. The standard deviations for error at a given value of x
    avg_stddev_data: The constant standard deviation (experimental error) that would justify the given fit.
    '''
    Y=f(Xdata,popt,dict_data)
    listdY=[]
    for i in xrange(len(popt)):
        p=popt[i]
        dp=abs(p)/1e6 + 1e-20
        popt[i]+= dp
        Yi=f(Xdata,popt,dict_data)
        dY=(Yi-Y)/dp
        listdY.append(dY)
        popt[i]-=dp
    listdY=scipy.array(listdY)
    #listdY is an array with N rows and M columns N=len(popt),M=len(Xdata[0])
    #pcov is an array with N rows and N columns
    left=scipy.dot(listdY.T,pcov)
    #Left is an array with M rows and N columns
    right=scipy.dot(left,listdY)
    #right is an M x M matrix. The diagonals of this are what we need.
    sigma2y=right.diagonal()
    #sigmay is the standard error of fit and is a function of x
    mean_sigma2y=scipy.mean(right.diagonal())
    M=Xdata.shape[1]
    N=len(popt)
    avg_stddev_data=scipy.sqrt(M*mean_sigma2y/N)
    #This is because if experimental error is constant at sig_dat, then mean_sigma2y=N/M*sig_dat**2
    sigmay=scipy.sqrt(sigma2y)
    return sigmay,avg_stddev_data
    popt=numpy.array([0.0]*11)
    
    
def fitdata(f,Xdata,Ydata,Errdata,pguess,dict_data,ax=False,ax2=False):
    '''
    fitdata(f,Xdata,Ydata,Errdata,pguess,dict_data)
        f=function f(X,p,dict_data)
        Xdata=array like object (k,M) shaped array for data with k predictors
        e.g. if X = (x1,x2,x3) then X=(X1,X2,X3) where X1 is a vector of x1 etc
        Ydata=array like object of length M
        Errdata=array like object of length M: error estimate of ydata.
        pguess=array like object of length N(vector of guess of parameters)
        dict_data= dictionary containing other data necessary for f
    Returns:
        popt=vector of length N of the optimized parameters
        pcov=covariance matrix of the fit
        perr=vector of length N of the std-dev of the optimized parameters
        p95=half width of the 95% confidence interval for each parameter i.e. popt-p95 and popt+p95
        p_p=vector of length N of the p-value for the parameters being zero
                (if p<0.05,null hypothesis rejected and parameter is non-zero)
        chisquared=(chisquared,chisquared_red,degfreedom,p)
        chisquared=chisquared value for the fit:sum of squared of weighted residuals
        chisquared_red=chisquared/degfreedom. Value should be approx. 1 for a good fit.
        degfreedom=M-N the degrees of freedom of the fitting
        chisquare=(p,chisquared,chisquared_red,degfreedom)
                p=Probability of finding a chisquared value at least as extreme as the one shown
                  purely by random chance(should be high for good fit)
                chisquared=chisquared value for the fit: sum of squares of weighted residuals
                chisquared_red=chisquared/degfreedom. Value should be  approx. 1 for a good fit.
                degfreedom=M-N the  degrees of freedom of the fitting
                R2=correlation coefficient or proportion of explained variance
                R2_adj=adjusted R2 taking into account number of predictors
        resanal=(p,w,mean,stddev) Analysis of residuals
                p=Probability of finding a w at least as extreme as the one observed (should be high for good fit)
                w=Shapiro-wilk test criterion
                mean=mean of residuals
                p_res=probability that the mean value obtained is different form zero merely by chance
                (should be low for good fit)
                The mean must be within 1 stddev of zro for highly significant fitting.
                F=F-statistic for the fit MSM/MSE
                Null hypothesis is that there is NO Difference between the twpo variances.
                p_F=probability that this value of F can arise by chance alone.
                   p_F<0.05 to reject null hypothesis and prove that the fit is good
                dw=Durbin_Watson statistic (value between 0 and 4).
                   2=no-autocorrelation.  0=+ve autocorrelation, 4 = -ve autocorrelation.
    '''
    
    def error(p, Xdata, Ydata, Errdata, dict_data):
        Y=f(Xdata, p,dict_data)
        residuals= (Y-Ydata)/Errdata
        return residuals
    
    res = scipy.optimize.leastsq(error, pguess, args=(Xdata, Ydata, Errdata, dict_data), full_output=1)
    (popt, pcov, infodict, errmsg, ier) = res
    perr=scipy.sqrt(scipy.diag(pcov))
    
    M=len(Ydata)
    N=len(popt)
    #Residuals
    Y=f(Xdata,popt,dict_data)
    residuals=(Y-Ydata)/Errdata
    meanY=scipy.mean(Ydata)
    squares=(Y-meanY)/Errdata
    squaresT=(Ydata-meanY)/Errdata
    
    SSM=sum(squares**2) #Corrected Sum of Squares
    SSE=sum(residuals**2) #Sum of Squares of Errors
    SST=sum(squaresT**2) #Total corrected sum of squares
    
    DFM=N-1 #Degrees of freedom for model
    DFE=M-N #Degrees of freedom for error
    DFT=M-1 #Degrees of freedom total
    
    MSM=SSM/DFM #Mean squares for model (explained variance)
    MSE=SSE/DFE #Mean squares for Error (should be small wrt MSM) Unexplained Variance
    MST=SST/DFT #Mean squares for total
    
    R2=SSM/SST #proportion of explained variance
    R2_adj=1-(1-R2)*(M-1)/(M-N-1) #Adjusted R2
    
    #t-test to see if parameters are different from zero
    t_stat=popt/perr #t-statistic for popt different from zero'
    t_stat=t_stat.real
    p_p=1.0-scipy.stats.t.cdf(t_stat,DFE) #should be low for good fit.
    z=scipy.stats.t(M-N).ppf(0.95)
    p95=perr*z
    #Chisquared Analysis on Residuals
    chisquared=sum(residuals**2)
    degfreedom=M-N
    chisquared_red=chisquared/degfreedom
    p_chi2=1.0-scipy.stats.chi2.cdf(chisquared,degfreedom)
    stderr_reg=scipy.sqrt(chisquared_red)
    chisquare=(p_chi2,chisquared, chisquared_red, degfreedom,R2,R2_adj)
    
    #Analysis of Residuals
    w,p_shapiro=scipy.stats.shapiro(residuals)
    mean_res=scipy.mean(residuals)
    stddev_res=scipy.sqrt(scipy.var(residuals))
    t_res=mean_res/stddev_res #t-statistic to test that mean_res is zero.
    p_res=1.0-scipy.stats.t.cdf(t_res,M-1)
        #if p_res <0.05, null hypothesis rejected and mean is non-zero.
        #Should be high for good fit.
    #F-test on residuals
    F=MSM/MSE #explained variance/unexplained . Should be large
    p_F=1.0-scipy.stats.f.cdf(F,DFM,DFE)
        #if p_F <0.05n, null-hypothesis is rejected
        #i.e. R^2>0 and at least one of the fitting parameters >0.
    dw=stools.durbin_watson(residuals)
    resanal=(p_shapiro,w,mean_res,F,p_F,dw)
    
    
    
    if ax:
        formataxis(ax)
        ax.plot(Ydata,Y,'ro')
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
        titletext='Parity plot for fit.\n'
        titletext+=r'$r^2$=%5.2f, $r^2_{adj}$=%5.2f, '
        titletext +='$\sigma_{exp}$=%5.2f,$\chi^2_{\nu}$=%5.2f, $p_{\chi^2}$=%5.2f, '
        titletext +='$\sigma_{err}^{reg}$=%5.2f'
        
        ax.title.set_text(titletext%(R2,R2_adj, avg_stddev_data, chisquared_red,  p_chi2, stderr_reg))
        ax.figure.canvas.draw()
        
    if ax2:#Test for homoscedasticity
        formataxis(ax2)
        ax2.plot(Y,residuals,'ro')
        
        ax2.xaxis.label.set_text('Fitted Data')
        ax2.yaxis.label.set_text('Residuals')
        
        titletext='Analysis of Residuals\n'
        titletext+=r'mean=%5.2f, $p_{res}$=%5.2f, $p_{shapiro}$=%5.2f, $Durbin-Watson$=%2.1f'
        titletext+='\nF=%5.2f, $p_F$=%3.2e'
        ax2.title.set_text(titletext%(mean_res, p_res, p_shapiro, dw , F, p_F))
        
        ax2.figure.canvas.draw()
        
    return popt,pcov,perr,p95,p_p,chisquare,resanal
    
    
def import_data(xlfile,sheetname):
    df=pd.read_excel(xlfile,sheetname=sheetname)
    return df
    
def prepare_data(df,Criterion,Predictors,Error=False):
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
    return Xdata,Y,Errdata
    
if  __name__ =="__main__":
    fig=plt.figure()
    ax=fig.add_subplot(111)
    fig.show()
    
    fig2=plt.figure()
    ax2=fig2.add_subplot(111)
    fig2.show()
    
# Make arbitrary function of three variables
    def f(X,p,dict_data):
#        a=dict_data['a']
#        b=dict_data['b']
        (x,y)=X
        Y=p[0]+p[1]*x**2+p[2]*y
        return Y
    
    
    #Get data from excel file using pandas
    dict_data=[]
    df=import_data('SynthData5.xlsx','Data')
    Xdata,Ydata,Errdata = prepare_data(df,'Ydata',('x','y'),Error= 'err')
    
    #Initial Guess
    N=3
    pguess=N*[0.0]
    
    resanal=fitdata(f, Xdata, Ydata, Errdata, pguess,dict_data, ax=ax, ax2=ax2)
    print resanal  
            