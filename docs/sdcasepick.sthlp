{smcl}
{ul on}{bf:Title}{ul off}

   {bf:sdcasepick}  {c -}{c -}  Second difference calculation based on units with similar baseline probabilities

   
{ul on}{bf:Syntax}{ul off}

   {bf:sdcasepick }{it:treatment} [{it:if}] [{it:in}]{bf:,} {bf:{ul on}g{ul off}roupvar(}{it:varname}{bf:)} [{bf:{ul on}prv{ul off}ar(}{it:varname}{bf:)} 
   {bf:{ul on}pru{ul off}b(}#{bf:)} {bf:{ul on}prl{ul off}b(}#{bf:)} {bf:{ul on}b{ul off}asevalue(}#{bf:)} {bf:{ul on}re{ul off}fvalue(}#{bf:)} {bf:{ul on}ro{ul off}und(}#{bf:)}  {bf:{ul on}e{ul off}xpression(}{it:exp}{bf:)} 
   {bf:{ul on}ou{ul off}tcome(}#{bf:)} {bf:{ul on}l{ul off}evel(}#{bf:)} {bf:{ul on}s{ul off}econddifference {ul on}v{ul off}erbose} {bf:{ul on}d{ul off}ropprvar} {bf:{ul on}i{ul off}vlist(}{it:varlist}{bf:)} 
   {bf:{ul on}ov{ul off}ergroup {ul on}oo{ul off}s}]


{ul on}{bf:Description}{ul off}

   After the user estimates a logit (or other binary choice) 
   regression that includes a binary group indicator
   (or a categorical variable if two of its categories are the groups of interest)
   {bf:sdcasepick} allows the user to specify any covariate of 
   interest z included in the previously estimated regression and 
   a range of baseline probabilities; within that range, it finds 
   the pair of observations (one in each group) in the estimation 
   sample that are closest in terms of baseline probability,
   defined as Pr(Y=1|{bf:x},z=baseline value).
   
   It can also be used after estimating a mutinomial logit, in which case 
   the relevant quantity is Pr(Y=m|{bf:x},z=baseline value).
      
   Optionally, it then calculates first and second differences 
   (and associated standard errors) by setting the other covariates 
   at the values in this closest pair of observations.


{ul on}{bf:Required inputs}{ul off}

   {bf:treatment} : the covariate of interest, z.
   
   {bf:groupvar(}{it:varname}{bf:)} : the binary (0/1) group indicator.  This option must be 
   specified.  This variable need not be in the previously estimated model; 
   see {bf:Examples} for illustration.


{ul on}{bf:Optional inputs}{ul off}

   {bf:prvar(}{it:varname}{bf:)} : the variable containing the predicted
   probabilities of interest for each in-sample observation, typically
   Pr(Y=1|{bf:x},z). (After mlogit: Pr(Y=m|{bf:x},z).) If not specified, 
   a variable will be created.
      
   {bf:prub(}{it:#}{bf:)} : the upper bound of permissible baseline probabilities.
   If not specified, the default value is 1.
   
   {bf:prlb(}{it:#}{bf:)} : the lower bound of permissible baseline probabilities.  
   If not specified, the default value is 0.
   
   {bf:basevalue(}{it:#}{bf:)} : the value at which to set the covariate of interest
   z when calculating baseline probability. If not specified, the default value is 0. 
   
   {bf:refvalue(}{it:#}{bf:)} : the value at which to set the covariate of interest
   z when calculating the change from z=basevalue, i.e., the first difference. 
   If not specified, the default value is 1.   
   
   {bf:round(}{it:#}{bf:)} : the units in which to round z for the purpose of 
   identifying the pair of observations that are closest in baseline probability.  
   This is useful for continuous z, when there are few or no observations with {it:exactly} 
   the same value of z across groups.  Note that the values of {bf:basevalue} and {bf:refvalue}  
   are not rounded {it:when calculating second differences}, even if {bf:round} is specified.
   Entering {bf:round(.1)} rounds 1.234 to 1.2, {bf:round(.01)} to 1.23. See Stata's {bf:f_round}.
   If not specified, z is not rounded.
   
   {bf:expression(}{it:exp}{bf:)} : specify if calculating probabilities based on anything
   OTHER than the default expression used by {bf:margins} after the estimated regression.
   Rarely needed, except with firthlogit. See the option {bf:expression} with {bf:margins} for specifics.
   See also the option {bf:outcome}, which is mutually exclusive, and used after mlogit.
   
   {bf:outcome(}{it:#}{bf:)} : the value of m for which to calculate Pr(Y=m|{bf:x}) after a multinomial 
   logit.  I.e., this should be a (non-baseline) value of the dependent variable.
   This option must be specified after mlogit and is not allowed otherwise.
   
   {bf:level(}{it:#}{bf:)} : the confidence level for the confidence intervals on the 
   first and second differences.  If not specified, the level set in Stata 
   (by default, 95) is used.
   
   {bf:seconddifference} : calculate second differences. 
   
   {bf:verbose} : show output from margins.
   
   {bf:dropprvar} : drop the program-created variable storing predicted probabilities. 
   
   {bf:ivlist(}{it:varlist}{bf:)} : a list of covariates in the previously estimated regression,   
   including {bf:groupvar}. Since {bf:sdcasepick} can determine this list automatically,  
   this option is typically not needed.

   {bf:oos} : not currently documented. do not use.
   
   {bf:overgroup} : not currently documented.  do not use.
   

{ul on}{bf:Examples}{ul off}

   {bf:logit dv gv##(z x1 x2)}
   {bf:sdcasepick z, groupvar(gv) prlb(.45) prub(.55) seconddifference}

   {bf:logit dv gv##(z x1 x2)}
   {bf:predict pv}
   {bf:sdcasepick z, groupvar(gv) prlb(.45) prub(.55) prvar(pv) seconddifference}
   
   Using option round with continuous treatment:
   {bf:logit dv gv##(c.z c.x1 c.x2)}
   {bf:sdcasepick z, groupvar(gv) basevalue(1.2) refvalue(1.6) round(.1) seconddifference}
   
   Comparing effects across two levels of a categorical group variable:
   {bf:logit dv i.cgv##(z x1 x2)}
   {bf:gen gv = 0 if cgv == 3}
   {bf:replace gv = 1 if cgv == 6}
   {bf:sdcasepick z, groupvar(gv)  prlb(.45) prub(.55) seconddifference}
   
   Specification after mlogit:
   {bf:mlogit dv gv##(z x1 x2)}
   {bf:predict pv, outcome(3)}
   {bf:sdcasepick z, groupvar(gv) outcome(3) prvar(pv) prlb(.10) prub(.20) seconddifference}

   
{ul on}{bf:Reference}{ul off}

   Gregory A. Caldeira and Daniel Lempert. (2018). 
   "Selection of Cases for Discussion: the U.S. Supreme Court, OT 1939, 1968, and 1982." 
   Paper presented at 2018 APSA Annual Meeting. Forthcoming, {it:Journal of Law & Courts}.
   

   
{ul on}{bf:Notes}{ul off}   
	
   {bf:sdcasepick} currently requires Stata 15 or newer; possible to change to earlier version
   at/near top of .ado file if desired (but output from older versions has not been tested).
   
   In addition to logit and mlogit, {bf:sdcasepick} can also be used after probit, firthlogit, and
   scobit, though these models (esp. scobit) have been less thoroughly tested.  Firthlogit requires the
   correct expression for Pr(Y=1) to be given in the option expression, since that model's
   margins default is specified incorrectly.  Support for ologit may be forthcoming at some point,
   but the advantages over the LPM are probably minimal.
      
   {bf:sdcasepick} incorporates a very slightly modified version of Long and Freese's mlincom
   command and subroutines thereof.  See J. Scott Long and Jeremy Freese. (2014). 
   Regression Models for Categorical Dependent Variables Using Stata, Third Edition. 
   College Station, TX: Stata Press.  

   
{ul on}{bf:Author}{ul off}

   Daniel Lempert
   lemperds@potsdam.edu / dalempert@gmail.com
