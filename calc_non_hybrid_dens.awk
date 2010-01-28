BEGIN{
bin=0
sum=0
}

{

if ($1 ~ /^[0-9]/){
	d=fabs($1-adressc)
    if (d<adressw || d > (adressw+adressh)){
#if ($1 <3.046-delta || ($1 >4.046+delta && $1 <8.046-delta) || $1 > 9.046+delta){

                sum+=$2
                bin+=1
          
        }
}


}
END{
#printf "Average density: %f from %d\n", sum/bin, bin
print sum/bin
}

function fabs ( x ) { return (x >= 0) ? x : -x }
