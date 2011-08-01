#! /usr/bin/awk -f
BEGIN{
	#Version check, use gawk only !
	"awk --version" | getline version;
	if ( version !~ /GNU Awk/)
	{
		print "#WARING: This is not gawk (GNU awk). Could cause problem!";
	}
	print "#This is avg_bl 0.1.1 (16.12.09)";
	if (ARGC>2)
	{
		print "#WARING: Average over multiple files!!";
	}
	if (name=="")
	{
		if (ARGV[1]=="")
		{
			print "#FILE: STDIN";
		}
		else
		{
			print "#FILE:",ARGV[1];
		}
	}
	else
	{
		print "#FILE:",name,"(set by user)";
	}	
	if (col=="")
	{
		col=1;
		print "#Setting column to",col,"(use option \"-v col=X\" to change)";
	}
	else
	{
		print "#Using column",col,"(set by user)";
	}
	if (N_B=="")
	{
		N_B=16;
		print "#Setting Number of blocks to",N_B,"(use option \"-v N_B=X\" to change)";
	}
	else
	{
		print "#Using ",N_B," blocks ";
	}
	c=0;
	warn=0;
}

/^[#@]/{
	next;
}

{
	if (NF>=col)
	{
		x[c]=$col;
		c++;
	}
	else
	{
		if (warn == 0)
		{
			print "# Not enough data in line",NR;
			warn++;
		}
	}	
	
}

END{    {if (c <= N_B && c>1); N_B=c;
         if (c <= 1); die "Not enough data points"}
	L_B=int(c/N_B);
	sum=0;
	sum2=0;
	for (b=0;b<N_B;b++){
		for (i=0;i<L_B;i++){
			nr=int(b*L_B+i)
			sum+=x[nr];
			sum2+=x[nr]*x[nr];
			sum_b[b]+=x[nr];
		}
	}	
	sum=sum/(L_B*N_B);
	sum2=sum2/(L_B*N_B);
	se_n=(sum2-sum*sum)/(L_B*N_B);
	se_n=sqrt(se_n);
	se_b=0;
	for (b=0;b<N_B;b++){
		sum_b[b]=sum_b[b]/L_B;
		#print b,sum,sum_b[b],se;
		se_b+=(sum_b[b]-sum)*(sum_b[b]-sum);
	}
#	print se;
	se_b=se_b/(N_B-1)/N_B;
#	print se;
	se_b=sqrt(se_b);
	print "#Output format: avg se naiv_se values";
	print sum,se_b,se_n,L_B*N_B,c;
}
