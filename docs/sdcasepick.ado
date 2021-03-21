//version 0.6--12/25/18.--dl
//compared to 0.5:
	//added outcome option and other modifications to work with mlogit
	//fixed expression option to also apply to program-created prvar
	//removed indeplist requirement and related improvments, incl: 
		//making option ivlist superfluous
		//works (a bit theoretically) with scobit
//compared to 0.4:
	//adds level option, incorporates dllimcom/dl_rm_lincom_stats instead of requiring spost13
	//fixed major bug with round option
	//added random integer to name of var generated if prvar not specified
//compared to 0.3:
	//relaxes requirement that groupvar be IV in regression.
	//fixes parse error if option specified in regression w/o ivlist
//compared to 0.2:
	//adds option round
//compared to 0.1:
	//basevalue fixed so that equality holds with float precision 
//compared to 0.0: 
	//added noi to display Note for prvar.
	//added `basevalue' `refvalue' suffix to prvar to avoid problem w/ repeated calls on v fast computers
	//added dropprvar option.
//to do: 
	//!ensure robustness of effect calculation to order interaction term entered into in regression, when *only* treatvar and groupvar are interacted...
	//...noting that mlincom 4-3 gives different outputs after logit y x##z versus logit y z##x.
	//-allow gv to be categorical but specify two values to use? 
	//-add ci option (probably have to do thru contrast rather than lincom/dlmlincom).
//note: see msdcasepick.ado

cap program drop sdcasepick

program define sdcasepick, rclass
version 15

//could change to version 13 at some point if robust

syntax varname [if] [in], Groupvar(varname numeric fv)  [PRVar(varname numeric) Ivlist(varlist numeric fv) PRUb(real 1) ///
PRLb(real 0) Basevalue(real 0) REfvalue (real 1)  ROund(real 0) ///
Seconddifference OVergroup OOs  Verbose Expression(passthru) Level(passthru) OUtcome(passthru) Dropprvar]

//add "absolute" option if care about absolute effect size difference?... 
//...lincom does not work w/ abs() so would need to do some kind of -1*1 if 1<0 cond etc.
//note that base(0) other(1) gives slightly different than base(1) other(0)b/c different baseline observations 
//selected (by def'n of "base").
//probably delete overgroup option.
//think about if/in marksample... oos option?--check why seems to give same margins?
//some kind of error if no interaction included in regression?


**factor variable notes:
//note that fvrevar does not give error if var is not factor var,
//so fvrevar(x1) gives "x1".

quietly{

tempvar tp
tempvar n
tempname activereg
*tempvar gv

marksample touse
*cap drop tu
*gen tu=`touse'

if "`oos'"!="oos"{
local es= "& e(sample)"
}

if "`oos'"=="oos"{
local noes= "noesample"
}


if "`seconddifference'"=="seconddifference"{
est store `activereg' //to restore before program exits.
}


*set up multinomial logit specific options:
if ( "`outcome'" =="" & e(cmd)=="mlogit" )| ("`outcome'" !="" & e(cmd)!="mlogit"){ //be careful about 2nd condition if ologit etc added.
	noi di "{error: specify option outcome if and only if estimation command in memory is mlogit}"
	exit 198
}

if "`expression'" !="" & e(cmd)=="mlogit" {
	noi di "{error: option expression not allowed with mlogit; see option outcome}"
	exit 198
}

if e(cmd)=="mlogit"{
	local mnlpredict = "predict(`outcome')" //used in margins
}




//note local `varlist' gives varname, i.e., var whose effect user wants.

*di wordcount("`ivlist'")

gen `n'=_n if `touse' `es'

local ec `e(cmdline)'
local ec: subinstr local ec  "##" "", all
local ec: subinstr local ec  "#" "", all
local ec: subinstr local ec  "(" " ", all
local ec: subinstr local ec  ")" " ", all
local ec: subinstr local ec  "," " ", all

tokenize `ec',parse(" ")



if "`prvar'"==""{
	local time=subinstr(c(current_time),":","",.)
	local date=subinstr(c(current_date)," ","",.)
	local rint=runiformint(0,10000)
	noi di "Note: option prvar not specified so default predicted values created"
	if "`dropprvar'"==""{
		noi di "and stored in varname pr`time'`date'_`rint'"
	}
	qui predict pr`time'`date'_`rint' if `touse' `es' , `expression' `outcome'
	local prvar "pr`time'`date'_`rint'"
}

if "`groupvar'"!=""{
	fvrevar `groupvar', list 
	local gv  "`r(varlist)'"
	egen `tp'=tag(`prvar' `gv') if `touse' `es'
}

if "`groupvar'"==""{
	fvrevar `3', list //prob. could fvrevar ec upfront.
	local gv  "`r(varlist)'"
	egen `tp'=tag(`prvar' `gv') if `touse' `es'
	*egen `tp'=tag(`prvar' `3' ) //**add `varlist' in tag()?? prob not needed tho prob no harm?
	//where `3' is group var, assuming user
	//has entered regression as specified.
	//note that `prvar' does what we want here given contents of <if "`prvar'"==""">
	//loop immediately above.
}

*give error if `gv' not 0/1:
*noi tab `touse' if e(sample)
*noi tab `gv' if `touse' `es'
*noi tab `gv' if `touse'
*noi tab `gv' if e(sample)
levelsof `gv' if `touse' `es', local(ggg)
if "`ggg'"!="0 1"{
	noi di "{error:Group Indicator (`gv') must be 0/1}"
	exit 198
}	




levelsof `n' if `tp'==1 & `gv'==0 & round(`varlist',`round')==round(float(`basevalue'),`round') ///
& `prvar' < `prub' & `prvar' > `prlb' ,local(g0obs) //already touse/es in def'n of `n'

levelsof `n' if `tp'==1 & `gv'==1 & round(`varlist',`round')==round(float(`basevalue'),`round') ///
& `prvar' < `prub' & `prvar' > `prlb' ,local(g1obs)

if "`g0obs'"==""{ 
	noi di "{error:Bound constraint could not be met for `gv'==0. Try different values of prub or prlb, or see option round.}"
	exit 198
}

if "`g1obs'"==""{ 
	noi di "{error:Bound constraint could not be met for `gv'==1. Try different values of prub or prlb, or see option round.}"
	exit 198
}




local min_dist=1
foreach g0o of local g0obs{
	foreach g1o of local g1obs{
		local d=abs(`prvar'[`g0o']-`prvar'[`g1o'])
		if `d'<`min_dist'{
			local min_g0o=`g0o' 
			local min_g1o=`g1o' 
			local min_dist=`d'
		}
	}
} //already touse/es in def'n of `n' --> g0/1obs






if "`seconddifference'"!="seconddifference"{ 
noi di  _newline(1) "{result:{ul on}CLOSEST PAIR OF OBSERVATIONS{ul off}}" _newline(2) ///
"{res:Group==0: `min_g0o'}"  _newline(1) ///
"{res:Group==1: `min_g1o'}"  _newline(1) ///
"{res:Group Variable: `gv'}" _newline(1) ///
"{res:Distance (Prediction) Variable: `prvar'}" _newline(1) ///
"{res:Distance Between Obss: `min_dist'}"
}


if "`seconddifference'"=="seconddifference"{

	if wordcount("`ivlist'")!=0{ //i.e., if user entered list of ivs.
		local ivlist: list uniq local ivlist
		local ivs
		foreach v of local ivlist{
			local v=regexr("`v'","[0-9]+\.","") //remove all level prefixes.
			*di "`v'"
			if regexm("`v'","#")==0{ //only non-interactions
				local ivs "`ivs' `v'" 
				*di "`ivs'"
			}
		}
		local ivs: list uniq local ivs //remove duplicates again once level prefixes removed.
		local ivs: list  local ivs - local varlist	//remove treatment	
		//nb: keep gv in list if in as iv!
			
		local m: word count `ivs' //returned from fvrevar
		local eq "=" //admittedly awkward
		forvalues i=1/`m'{
			*local g0x`i'=``i''[`min_g0o'] //i.e., set at observed value for _n=min_g0o
			*local g1x`i'=``i''[`min_g1o']
			local g0x`i'= `: word `i' of `ivs''[`min_g0o'] //i.e., set at observed value for _n=min_g0o
			local g1x`i'= `: word `i' of `ivs''[`min_g1o']
			local atvals0 `atvals0' `: word `i' of `ivs'' `eq' `g0x`i'' 
			local atvals1 `atvals1' `: word `i' of `ivs'' `eq' `g1x`i'' 
		}
	}
		
	if wordcount("`ivlist'")==0{ //i.e., if user did not enter list of ivs.
		//code to get list of rhs vars--addition/mod to Nick Cox's code
		//https://www.stata.com/statalist/archive/2006-03/msg00372.html
		mat B = e(b) 
		tokenize "`: colnames e(b)'" 
		local rhsvlist
		forval j = 1/`= colsof(B)' { 
			if B[1, `j'] != 0 & "``j''" != "_cons"  & "``j''" != "lnalpha" { 
				//lnalpha condition added for scobit
				local rhsvlist `rhsvlist' ``j''
			}
		} 
		di "`rhsvlist'"
		local rhsvlist: list uniq local rhsvlist
		local ivs

		
		foreach v of local rhsvlist{
			local v=regexr("`v'","[0-9]+\.","") //remove all level prefixes.
			*di "`v'"
			if regexm("`v'","#")==0{ //only non-interactions
				local ivs "`ivs' `v'" 
				*di "`ivs'"
			}
		}
		local ivs: list uniq local ivs //remove duplicates again once level prefixes removed.
		local ivs: list  local ivs - local varlist		
		//nb: keep gv in list if in as iv!

		
		local m: word count `ivs'
		local eq "=" //admittedly awkward
		forvalues i=1/`m'{
			local g0x`i'= `: word `i' of `ivs''[`min_g0o'] //i.e., set at observed value for _n=min_g0o
			local g1x`i'= `: word `i' of `ivs''[`min_g1o']
			local atvals0 `atvals0' `: word `i' of `ivs'' `eq' `g0x`i'' 
			local atvals1 `atvals1' `: word `i' of `ivs'' `eq' `g1x`i'' 
			*di "`atvals0'"
		}
	}
		
	if "`verbose'"=="verbose"{
		local noi noisily
	}
	
	if "`overgroup'" !="overgroup"{
		`noi' margins if `touse', at(`atvals0' `varlist'=(`basevalue' `refvalue'))  ///
		at( `atvals1' `varlist'=(`basevalue' `refvalue')) post `noes' `expression'  `mnlpredict'
		//note that `gv' not needed above b/c `atvals0/1' incorporates `gv' (or cat. covariate on which its based).
		`noi' qui dlmlincom 2-1, `level' all rowname(Effect:Group 0) estname(effect) clear 
		ret scalar g0_level=`r(level)'
		ret scalar g0_ub=`r(ub)'
		ret scalar g0_lb=`r(lb)'
		ret scalar g0_p=`r(p)'
		ret scalar g0_z=`r(z)'
		ret scalar g0_se=`r(se)'
		ret scalar g0_est=`r(est)'
		`noi' qui dlmlincom 4-3, `level' all rowname(Effect:Group 1) estname(effect) add 
		ret scalar g1_level=`r(level)'
		ret scalar g1_ub=`r(ub)'
		ret scalar g1_lb=`r(lb)'
		ret scalar g1_p=`r(p)'
		ret scalar g1_z=`r(z)'
		ret scalar g1_se=`r(se)'
		ret scalar g1_est=`r(est)'
		noi di  _newline(1) "{result:{ul on}EFFECT SIZES BY GROUP AND SECOND DIFFERENCE (`varlist') {ul off}}"
		*the following two ifs can be removed and can leave single line "noi lincom..." if
		*long fixes dlmlincom to allow level option.  
		*until then, require dldlmlincom if allow level option in sdcasepick
		noi dlmlincom (4-3)-(2-1), `level' all rowname(Effect:2nd Diff) estname(effect) add // detail
		ret scalar sd_level=`r(level)'
		ret scalar sd_ub=`r(ub)'
		ret scalar sd_lb=`r(lb)'
		ret scalar sd_p=`r(p)'
		ret scalar sd_z=`r(z)'
		ret scalar sd_se=`r(se)'
		ret scalar sd_est=`r(est)'
	}
	
	if "`overgroup'" =="overgroup"{ //if keep in, check robustness carefully.
		`noi' margins, at(`atvals0' `varlist'=(`basevalue' `refvalue'))  ///
		at( `atvals1' `varlist'=(`basevalue' `refvalue' )) over(`gv') ///
		post `noes' `expression'  `mnlpredict'
		qui dlmlincom 3-1, `level' all rowname(Effect:Group 0) clear
		qui dlmlincom 8-6, `level' all rowname(Effect:Group 1) add
		noi dlmlincom (8-6)-(3-1), `level' all rowname(Effect:2nd Diff)  ///
		title("Group Effects and Second Difference") add
	}
		
	noi di  _newline(1) "{result:{ul on}CLOSEST PAIR OF OBSERVATIONS{ul off}}" _newline(2) ///
	"{res:Group==0: `min_g0o'}"  _newline(1) ///
	"{res:Group==1: `min_g1o'}"  _newline(1) ///
	"{res:Group Variable: `gv'}" _newline(1) ///
	"{res:Distance (Prediction) Variable: `prvar'}" _newline(1) ///
	"{res:Rhs Vars: `ivs'}" _newline(1) ///
	"{res:Distance Between Obss: `min_dist'}"
	
	est restore `activereg'	
}	
	
if "`dropprvar'"=="dropprvar"{
cap drop pr`time'`date'`pvbv'`pvov' //cap in case prvar (nonsensically) also specified
}
	
capture ret local gv = "`gv'"
capture ret local prvar = "`prvar'"
capture ret scalar dist=`min_dist'
capture ret scalar g0_obs=`min_g0o'
capture ret scalar g1_obs=`min_g1o'
capture ret local g0_atvals="`atvals0'"
capture ret local g1_atvals="`atvals1'"
capture ret local ivs = "`ivs'"
capture ret local expression = "`expression'"
capture ret local treatment="`varlist'"
capture ret local basevalue=`basevalue'
capture ret local refvalue=`refvalue'

} //end quietly	
end




//**************************************************************************************************************
//----------------BELOW: SLIGHT MODIFICATIONS TO LONG/FREESE'S _rm_lincom_stats AND mlincom PROGRAMS------------

//dl--12/22/18: v. slight modification to long and freese's _rm_lincom_stats
//changes explained with "--dl:" prefix below  
// version 1.1.0 2016-09-14 | long freese | t-value for svy
// version 1.0.0 2014-02-14 | long freese | spost13 release

//  compute statistics from lincom returns

cap program drop dl_rm_lincom_stats

program define dl_rm_lincom_stats, rclass
version 11.2


    tempname estval seval levelval levelpval cifactorval ///
        zval pval lbval ubval df

    scalar `df' = r(df) // 2016-09-14
    scalar `estval'     = r(estimate) // lincom estimate
    scalar `seval'      = r(se) // lincom std err
    scalar `zval'       = `estval'/`seval' // zvalue

    if `df' < . { // 2016-09-14

        scalar `pval'       = 2*(1-t(`df',abs(`zval'))) // pvalue 2 tailed
        scalar `levelval'   = r(level) // --dl: modified--level for CI
        scalar `levelpval'  = ( 1 - (`levelval'/100) )/2 // pval for ci
        scalar `cifactorval' = abs(invt(`df',`levelpval')) // +- ci multiplier

    }
    else {

        scalar `pval'       = 2*(1-normal(abs(`zval'))) // pvalue 2 tailed
        scalar `levelval'   = r(level) // --dl: modified--level for CI
        scalar `levelpval'  = ( 1 - (`levelval'/100) )/2 // pval for ci
        scalar `cifactorval' = abs(invnormal(`levelpval')) // +- multiplier for ci

    }

    scalar `lbval'      = `estval' - (`cifactorval'*`seval') // lb
    scalar `ubval'      = `estval' + (`cifactorval'*`seval') // ub

    return scalar est   = `estval'
    return scalar se    = `seval'
    return scalar z     = `zval'
    return scalar p     = `pval'
    return scalar lb    = `lbval'
    return scalar ub    = `ubval'
    return scalar level = `levelval'

end



//dl--12/22/18: slight modification to long and freese's mlincom
//changes explained with "--dl:" prefix below  
* version 1.0.3 2015-01-09 | long freese | tweak row margin
* version 1.0.2 2014-10-01 | long freese | always allow non margins
* version 1.0.1 2014-09-28 | long freese | force if not margins
* version 1.0.0 2014-02-18 | long freese | spost13 release

//  generate lincom expression and sent it to lincom

*   DO: trap errors from lincom

cap program drop dlmlincom

program define dlmlincom

    version 11.2
    tempname newmat newmatall

    gettoken token 0 : 0, parse(",= ")
    while `"`token'"'!="" & `"`token'"'!="," {
        local expression `"`expression'`token'"'
        gettoken token 0 : 0, parse(",= ")
    }
    if `"`0'"'!="" {
        local 0 `",`0'"'
        syntax [, ///
        /// results to display
        /// force
            STATS(string asis) STATistics(string asis) ///
            ALLstats /// all stats (p, z, se, level)
            Details /// show lincom output
            NOTABle /// only show lincom
			LEVEL(passthru) /// --dl: pass level to lincom
        /// save results to matrix
            clear /// clear matrix before saving
            add save /// add results to matrix
        /// label matrix
            ROWName(string) label(string) /// row name for current estimate
            ROWEQnm(string) /// row eq name for current estimate
            ESTName(string asis) /// allow override margin name for estimate
        /// displaying matrix
            DECimals(integer 3) /// Digits when list matrix
            WIDth(integer 8) /// Column width when listing matrix
            title(string) /// Title when listing
            TWIDth(integer 0) /// 2015-01-09
        ]
    }
    else {
        local decimals 3
        local width 8
    }

    if ("`twidth'"=="0") local twopt ""
    else local twopt "twidth(`twidth')" // 2015-01-09

    if ("`label'"=="") local label `"`rowname'"' // synonyms
    if ("`save'"=="save") local add "add" // synonyms
    if ("`details'"=="details") local quietly ""
    else local quietly "quietly"
    * if no table, show lincom output
    if ("`notable'"=="notable") local quietly ""

    local matrix _mlincom
    capture confirm matrix `matrix'
    if (_rc == 0) local matrixexists = 1
    else local matrixexists = 0
    local matrixall _mlincom_allstats
    capture confirm matrix `matrixall'
    if (_rc == 0) local matrixallexists = 1
    else local matrixallexists = 0

    if "`expression'"=="" {
        if `matrixexists' == 1 {
            if "`clear'"=="clear" {
                capture matrix drop `matrix'
                capture matrix drop `matrixall'
                exit
            }
            * if no expression, list table
            matlist `matrix', format(%`width'.`decimals'f) title("`title'") ///
                `twopt'
            exit
        }
        else if ("`clear'"=="clear") exit
        else {
            display as error "you must specify the lincom expression"
            exit
        }
    }

    capture confirm matrix e(b)
    if _rc>0 {
        display as error ///
"mlincom requires e(b) in memory; with margins or mtable, use the post option"
        exit
    }

/*  dropped in 1.0.2 to allow mlincom to work with all estimation commands
    if "`e(cmd)'"!="margins" & "`force'"=="" {
        display as error ///
            "mlincom must be run immediately following margins, post"
        exit
    }
    capture confirm matrix e(b)
    if _rc>0 {
        display as error ///
            "mlincom requires margins or mtable with the post option"
        exit
    }
*/

//  decode expression

    local lc "`expression'"
    foreach c in ( ) + - {
        local lc = subinstr("`lc'","`c'"," `c' ",.)
    }

    * remove o. names from column names
    tempname b
    matrix `b' = e(b)
    local orignms : colfullnames `b'
    local noomitnms "" // without o. names
    foreach var of local orignms {
        _ms_parse_parts `var'
        if (!`r(omit)') local noomitnms `noomitnms' `var'
    }
    local lcstr ""
    foreach e in `lc' {
        if ("`e'"=="(") local lcstr "`lcstr'("
        else if ("`e'"==")") local lcstr "`lcstr')"
        else if ("`e'"=="+") local lcstr "`lcstr'+"
        else if ("`e'"=="-") local lcstr "`lcstr'-"
        else {
            local i = int(real("`e'"))
            if `i' == . {
                display as error "unexpected missing value found"
                exit
            }
            else {
                local bnm : word `i' of `noomitnms'
                local lcstr "`lcstr'_b[`bnm']"
            }
        }
    }

//  run lincom and compute stats

    `quietly' lincom `lcstr' , `level' //--dl: added "`level'"
    dl_rm_lincom_stats //--dl: modified version of _rm_lincom_stats
    scalar estimate = r(est)
    scalar se = r(se)
    scalar zvalue = r(z)
    scalar pvalue = r(p)
    scalar ll = r(lb)
    scalar ul = r(ub)
    scalar lvl = r(level)

//  add results to matrix

    * if not adding, clear matrix
    if "`clear'" == "clear" | "`add'"=="" {
        capture matrix drop `matrix'
        capture matrix drop `matrixall'
        local matrixexists = 0
    }

    local estimatenm "lincom"
    if ("`estname'"!="") local estimatenm "`estname'"

    * what stats go in table
    if ("`stats'"!="" & "`statistics'"=="") local statistics "`stats'"
    local statlist "estimate pvalue ll ul lvl" // --dl: default statistics--added lvl
    local statsall "estimate se zvalue pvalue ll ul lvl" //--dl: added lvl
    if ("`allstats'"=="allstats") local statlist "`statsall'"
    if ("`statistics'"=="noci") local statistics "est"
    if "`statistics'"=="all" {
        local statlist "`statsall'"
    }
    else if "`statistics'"!="" {
        local newstatlist ""
        foreach opt in `statistics' {
            _parse_stat `opt'
            local newopt
            local stat "`s(stat)'"
            local newstatlist "`newstatlist'`stat' "
            if "`s(isbad)'"=="1" {
                display as error ///
                    "invalid statistic specified: `opt'"
                exit
            }
        }
        local statlist "`newstatlist'"
    }
    if ("`s(isbad)'"=="1") exit
    local statlist : list uniq statlist

    * column names for matrix
    foreach stat in `statlist' {
        if "`stat'"=="estimate" { // use option estimatenm
            local `stat' "estimate"
            local colnms "`colnms'`estimatenm' "
        }
        else {
                local `stat' "`stat'"
            local colnms "`colnms'`stat' "
        }
    }
    local colnms = trim("`colnms'")

    * if add, make sure column names match
    if "`add'"=="add" {
        if `matrixexists' == 1 {
            local priorcolnms : colnames `matrix'
        }
        if `matrixexists'==1 & "`colnms'"!="`priorcolnms'" {
            display as error ///
                "statistics in matrix `matrix' do not match those being added"
            exit
        }
    }

    * get ALL statistics to be saved
    foreach s in `statsall' {
        local colnmsall "`colnmsall'`s' "
        matrix `newmatall' = nullmat(`newmatall') , `s'
    }
    * list of selected statistics s084
    foreach s in `statlist' {
        matrix `newmat' = nullmat(`newmat') , `s'
    }

    matrix colname `newmat' = `colnms'
    matrix colname `newmatall' = `colnmsall'
    if  "`roweqnm'" != "" {
        matrix roweq `newmat' = `roweqnm'
        matrix roweq `newmatall' = `roweqnm'
    }

    if "`label'"!="" { // nolabel
        matrix rowname `newmat' = `"`label'"'
        matrix rowname `newmatall' = `"`label'"'
    }
    else { // label
        local n = 1
        if (`matrixexists'==1) local n = rowsof(`matrix') + 1
        matrix rowname `newmat' = "`n'"
        matrix rowname `newmatall' = "`n'"
    }
    if `matrixexists'==1 { // it has been deleted if add not specified
        matrix `matrix' = `matrix' \ `newmat'
        matrix `matrixall' = `matrixall' \ `newmatall'
    }
    else {
        matrix `matrix' = `newmat'
        matrix `matrixall' = `newmatall'
    }

    if "`notable'"=="" {
        matlist `matrix', format(%`width'.`decimals'f) title("`title'") ///
            `twopt'
    }

end


//dl: _parse_stat is part of long and freese's mlincom
cap program drop _parse_stat
program define _parse_stat, sclass
    local isbad = 1
    local stat "`1'"
    local is = inlist("`stat'","e","es","est","esti","estim","estima")
    if `is'==1 {
        local stat "estimate"
        local isbad = 0
    }
    local is = inlist("`stat'","estimat","estimate","coef")
    if `is'==1 {
        local stat "estimate"
        local isbad = 0
    }
    local is = inlist("`stat'","s","se","stderr")
    if `is'==1 {
        local stat "se"
        local isbad = 0
    }
    local is = inlist("`stat'","p","pv","pva","pval","pvalu","pvalue")
    if `is'==1 {
        local stat "pvalue"
        local isbad = 0
    }
    local is = inlist("`stat'","z","zv","zva","zval","zvalu","zvalue")
    if `is'==1 {
        local stat "zvalue  "
        local isbad = 0
    }
    local is = inlist("`stat'","upper","ub","u","ul")
    if `is'==1 {
        local stat "ul"
        local isbad = 0
    }
    local is = inlist("`stat'","lower","lb","l","ll")
    if `is'==1 {
        local stat "ll"
        local isbad = 0
    }
    sreturn local stat "`stat'"
    sreturn local isbad "`isbad'"
end

