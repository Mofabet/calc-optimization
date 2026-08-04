#!/usr/bin/env bash
# ---------------------------------------------------------------------------
#  LMTO pipeline
#    1. ~/bin/lminit.soc   (retried until it finishes without an abort)
#    2. ~/bin/lmctl.soc
#    3. CTRL: IO VERBOS -> 50
#    4. ~/bin/lmctl.soc
#    5. CTRL: SCELL PLAT, vector i * FACTOR<i>   (skipped if all factors are 1)
#    6. ~/bin/lm47.run lmscell   (rewrites CTRL in place; skipped as well)
#    7. summary.dat: celldm, CELL_PARAMETERS, ATOMIC_POSITIONS, K_POINTS
#
#  Usage:  ./lmto_run.sh          (from the working directory)
#  Options via environment, e.g.:  MAXTRY=50 FACTOR3=2 SUMMARY=out.dat ./lmto_run.sh
# ---------------------------------------------------------------------------
set -uo pipefail

LMINIT=${LMINIT:-$HOME/bin/lminit.soc}
LMCTL=${LMCTL:-$HOME/bin/lmctl.soc}
LM47=${LM47:-$HOME/bin/lm47.run}
CTRL=${CTRL:-CTRL}            # lmscell overwrites this same file
MAXTRY=${MAXTRY:-30}          # max lminit.soc attempts
# supercell multipliers for the diagonal of SCELL PLAT.
# With all three equal to 1 the SCELL/lmscell steps are skipped completely.
FACTOR1=${FACTOR1:-1}
FACTOR2=${FACTOR2:-1}
FACTOR3=${FACTOR3:-1}
# add the last vertex of the SYML path to K_POINTS (1 = yes, 0 = no)
ADD_ENDPOINT=${ADD_ENDPOINT:-0}
SUMMARY=${SUMMARY:-summary.dat}
INITLOG=${INITLOG:-lminit.log}  # log of the successful lminit.soc run

say() { printf '=== %s\n' "$*"; }
die() { printf '\nERROR: %s\n' "$*" >&2; exit 1; }

# run a command quietly; on failure show the tail of its output and keep it
run_quiet() {                 # run_quiet <keep-on-error-file> <cmd> [args...]
    local errfile=$1; shift
    local tmp; tmp=$(mktemp ./.lmout.XXXXXX)
    "$@" > "$tmp" 2>&1
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        mv -f "$tmp" "$errfile"
        tail -n 20 "$errfile" >&2
        die "$1 failed (exit $rc), see $errfile"
    fi
    rm -f "$tmp"
}

# --------------------------- awk: patch SCELL ------------------------------
read -r -d '' SCELL_AWK << 'AWKEOF'
function isnum(x){ return (x ~ /^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([EeDdQq][+-]?[0-9]+)?$/) }
function process(   i,j,nf,fld,cnt,s,val,dirty,out,tok,new,row,fac){
  s=0; cnt=0
  for(i=1;i<=n;i++){
    nf=split(L[i],fld,/[ \t]+/); dirty=0
    for(j=1;j<=nf;j++){
      tok=fld[j]
      if(tok=="") continue
      if(s==2) continue
      if(s==0){
        # PLAT= may carry the first number glued to it
        if(tok ~ /^PLAT=/){
          s=1; cnt=0; val=substr(tok,6)
          if(val!="" && isnum(val)){                 # 1st number may be glued to PLAT=
            cnt=1; seen++
            if(F1!=1){
              new=sprintf("%.8f", val*F1)
              if(val+0 != 0) print val" > "new > CHG
              fld[j]="PLAT="new; dirty=1
            }
          }
        }
      } else {
        if(isnum(tok)){
          cnt++; seen++
          row=int((cnt-1)/3)+1                       # 1..3 = PLAT vector number
          fac=(row==1?F1:(row==2?F2:F3))
          if(fac!=1){
            new=sprintf("%.8f", tok*fac)
            if(tok+0 != 0) print tok" > "new > CHG   # report non-zero changes
            fld[j]=new; dirty=1
          }
          if(cnt==9) s=2
        } else s=2
      }
    }
    if(dirty){
      out=""
      for(j=1;j<=nf;j++){ if(fld[j]=="") continue; out=(out==""?fld[j]:out" "fld[j]) }
      if(L[i] ~ /^[ \t]/) out="        " out
      L[i]=out
    }
  }
  for(i=1;i<=n;i++) print L[i]
  n=0
}
BEGIN{ inb=0; n=0; seen=0 }
{
  if($0 ~ /^[^ \t]/){                    # a category starts in column 1
    if(inb){ process(); inb=0 }
    if($1 ~ /^SCELL/){ inb=1; n=0; L[++n]=$0; next }
  }
  if(inb){ L[++n]=$0; next }
  print
}
END{ if(inb) process(); if(seen!=9) exit 3 }
AWKEOF

# --------------------------- awk: summary ----------------------------------
read -r -d '' EXTRACT_AWK << 'AWKEOF'
function isnum(x){ return (x ~ /^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([EeDdQq][+-]?[0-9]+)?$/) }
{
  if($0 ~ /^[^ \t]/){ cat=$1; sub(/=.*/,"",cat) }   # category of this line
  if(cat=="STRUC") sbuf=sbuf" "$0
  else if(cat=="SITE") tbuf=tbuf" "$0               # ATOM= in CLASS is ignored
  else if(cat=="SYML") kbuf=kbuf" "$0
}
END{
  # ALAT and the 9 PLAT numbers from STRUC
  nf=split(sbuf,f,/[ \t]+/); alat=""; np=0
  for(i=1;i<=nf;i++){
    if(f[i]=="") continue
    if(f[i] ~ /^ALAT=/){ v=substr(f[i],6); if(v!="") alat=v; else { j=i+1; while(f[j]=="") j++; alat=f[j] } }
    if(f[i] ~ /^PLAT=/ && np==0){
      v=substr(f[i],6); if(v!="" && isnum(v)) plat[++np]=v
      for(j=i+1;j<=nf && np<9;j++){ if(f[j]=="") continue; if(!isnum(f[j])) break; plat[++np]=f[j] }
    }
  }
  # atom labels and positions from SITE (labels are taken as they are in CTRL)
  nf=split(tbuf,g,/[ \t]+/); na=0
  for(i=1;i<=nf;i++){
    if(g[i]=="") continue
    if(g[i] ~ /^ATOM=/){ na++; name[na]=substr(g[i],6); npos[na]=0 }
    if(g[i] ~ /^POS=/ && na>0){
      v=substr(g[i],5); if(v!="" && isnum(v)) pos[na,++npos[na]]=v
      for(j=i+1;j<=nf && npos[na]<3;j++){ if(g[j]=="") continue; if(!isnum(g[j])) break; pos[na,++npos[na]]=g[j] }
    }
  }
  if(alat=="" || np!=9 || na==0){
    print "cannot parse ALAT/PLAT/SITE (alat="alat" nplat="np" natom="na")" > "/dev/stderr"; exit 4
  }
  # high-symmetry path from SYML: NQ / Q1 / LAB1 of every segment
  nf=split(kbuf,h,/[ \t]+/); nk=0
  for(i=1;i<=nf;i++){
    if(h[i]=="") continue
    if(h[i] ~ /^NQ=/){ nk++; nq[nk]=substr(h[i],4); lab[nk]="" }
    else if(nk>0 && h[i] ~ /^Q1=/){
      v=substr(h[i],4); m=0; if(v!="" && isnum(v)) kp[nk,++m]=v
      for(j=i+1;j<=nf && m<3;j++){ if(h[j]=="") continue; if(!isnum(h[j])) break; kp[nk,++m]=h[j] }
    }
    else if(nk>0 && h[i] ~ /^LAB1=/) lab[nk]=substr(h[i],6)
    else if(nk>0 && h[i] ~ /^LAB2=/) endlab=substr(h[i],6)
    else if(nk>0 && h[i] ~ /^Q2=/){
      v=substr(h[i],4); m=0; if(v!="" && isnum(v)) endp[++m]=v
      for(j=i+1;j<=nf && m<3;j++){ if(h[j]=="") continue; if(!isnum(h[j])) break; endp[++m]=h[j] }
    }
  }

  printf "celldm = %s\n", alat
  print "CELL_PARAMETERS alat"
  for(r=0;r<3;r++) printf "%12.8f%12.8f%12.8f\n", plat[3*r+1], plat[3*r+2], plat[3*r+3]
  print "ATOMIC_POSITIONS alat"
  for(a=1;a<=na;a++){
    printf "%-6s", name[a]
    for(i=1;i<=npos[a];i++) printf "%12.8f", pos[a,i]
    printf "\n"
  }
  if(nk>0){
    print "K_POINTS"
    print nk + (ENDPT ? 1 : 0)
    for(k=1;k<=nk;k++)
      printf "%10.6f%10.6f%10.6f %4d !%s\n", kp[k,1], kp[k,2], kp[k,3], nq[k], lab[k]
    if(ENDPT) printf "%10.6f%10.6f%10.6f %4d !%s\n", endp[1], endp[2], endp[3], 1, endlab
  }
}
AWKEOF

# ---------------------------------------------------------------------------
# 1. lminit.soc -- repeat until it exits cleanly and produces CTRL
# ---------------------------------------------------------------------------
say "1. $LMINIT"
ok=0
for ((i=1;i<=MAXTRY;i++)); do
    tmp=$(mktemp ./.lminit.XXXXXX)
    "$LMINIT" > "$tmp" 2>&1
    rc=$?
    if [[ $rc -eq 0 ]] \
       && ! grep -qiE 'abort|core dumped|segmentation|forrtl' "$tmp" \
       && [[ -s $CTRL ]]; then
        mv -f "$tmp" "$INITLOG"     # keep only the successful run's log
        echo "   attempt $i: ok"
        ok=1; break
    fi
    rm -f "$tmp"
    echo "   attempt $i: abort > restart"
    sleep 1
done
[[ $ok -eq 1 ]] || die "lminit.soc did not succeed in $MAXTRY attempts"
[[ -s $CTRL ]] || die "$CTRL was not created"

# ---------------------------------------------------------------------------
# 2. lmctl.soc
# ---------------------------------------------------------------------------
say "2. $LMCTL"
run_quiet lmctl.err "$LMCTL"

# ---------------------------------------------------------------------------
# 3. IO VERBOS -> 50 (only inside the IO category)
# ---------------------------------------------------------------------------
say "3. IO VERBOS -> 50"
tmp=$(mktemp ./.lmctrl.XXXXXX)
awk '{ if($0 ~ /^[^ \t]/){ cat=$1 }
       if(cat=="IO") gsub(/VERBOS=[0-9]+/,"VERBOS=50")
       print }' "$CTRL" > "$tmp" && grep -q 'VERBOS=50' "$tmp" \
    || { rm -f "$tmp"; die "failed to set VERBOS=50 in $CTRL"; }
mv -f "$tmp" "$CTRL"

# ---------------------------------------------------------------------------
# 4. lmctl.soc again
# ---------------------------------------------------------------------------
say "4. $LMCTL"
run_quiet lmctl.err "$LMCTL"

# ---------------------------------------------------------------------------
# 5-6. supercell -- only if at least one factor differs from 1
# ---------------------------------------------------------------------------
STEP=5
if ! awk -v a="$FACTOR1" -v b="$FACTOR2" -v c="$FACTOR3" 'BEGIN{exit !(a==1&&b==1&&c==1)}'; then
STEP=7

say "5. SCELL: PLAT diagonal * $FACTOR1 $FACTOR2 $FACTOR3"
grep -q '^SCELL' "$CTRL" || die "no SCELL category in $CTRL"
tmp=$(mktemp ./.lmctrl.XXXXXX); chg=$(mktemp ./.lmchg.XXXXXX)
awk -v F1="$FACTOR1" -v F2="$FACTOR2" -v F3="$FACTOR3" -v CHG="$chg" "$SCELL_AWK" "$CTRL" > "$tmp"
if [[ $? -ne 0 ]]; then
    rm -f "$tmp" "$chg"
    die "could not patch PLAT in SCELL (9 numbers not found), $CTRL left untouched"
fi
mv -f "$tmp" "$CTRL"
while read -r line; do echo "   $line"; done < "$chg"
rm -f "$chg"

say "6. $LM47 lmscell"
run_quiet lmscell.err "$LM47" lmscell

fi

# ---------------------------------------------------------------------------
# 7. summary
# ---------------------------------------------------------------------------
say "$STEP. summary -> $SUMMARY"
awk -v ENDPT="$ADD_ENDPOINT" -f /dev/stdin "$CTRL" > "$SUMMARY" << EOF
$EXTRACT_AWK
EOF
[[ $? -eq 0 && -s $SUMMARY ]] || die "could not parse $CTRL"

say "Job done!"
