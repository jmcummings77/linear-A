#ifdef __cplusplus
#include <memory>
#endif
/* Shared protocol driver; native C++ calls its public CSRMatrix interface. */
static double sparse_token(void){double x;if(scanf("%lf",&x)!=1||!isfinite(x))fail("invalid sparse input token");return x;}
static size_t sparse_index(void){double x=sparse_token();if(x<0||x>=(double)SIZE_MAX||floor(x)!=x)fail("invalid CSR index");return (size_t)x;}
static int sparse_cli(int argc,char **argv){
    if(argc!=(argc>2&&!strcmp(argv[2],"gmres")?13:12))fail("expected sparse OP ROWS COLS NNZ ITERATIONS RTOL ATOL LIMIT JACOBI CAPTURE");
    const char *op=argv[2];if(strcmp(op,"spmv")&&strcmp(op,"dense")&&strcmp(op,"cg")&&strcmp(op,"gmres")&&strcmp(op,"ilu_setup")&&strcmp(op,"ilu_apply"))fail("invalid sparse operation");
    size_t rows=integer(argv[3]),cols=integer(argv[4]),nnz=integer(argv[5]),iterations=integer(argv[6]),limit=integer(argv[9]),jacobi=integer(argv[10]),capture=integer(argv[11]);
    double rtol=number(argv[7]),atol=number(argv[8]);if(jacobi>3||(strcmp(op,"gmres")&&jacobi>1)||capture>1||rows==SIZE_MAX||rows>SIZE_MAX/sizeof(size_t)-1||nnz>SIZE_MAX/sizeof(double))fail("invalid sparse dimensions/options");
    size_t *rp=(size_t*)calloc(rows+1,sizeof(size_t)),*ci=(size_t*)calloc(nnz?nnz:1,sizeof(size_t));double *v=(double*)calloc(nnz?nnz:1,sizeof(double));
    if(!rp||!ci||!v)fail("sparse allocation failed");
    for(size_t i=0;i<=rows;i++)rp[i]=sparse_index();for(size_t i=0;i<nnz;i++)ci[i]=sparse_index();for(size_t i=0;i<nnz;i++)v[i]=sparse_token();
    size_t count=(!strcmp(op,"cg")||!strcmp(op,"gmres"))?rows:cols;matrix b={0},dense={0};require(m_create(count,1,&b));read_values(&b);
    int trailing;do{trailing=getchar();}while(trailing!=EOF&&isspace((unsigned char)trailing));if(trailing!=EOF)fail("extra sparse input");
#ifdef __cplusplus
    linear_a::CSRMatrix sparse(rows,cols,std::vector<size_t>(rp,rp+rows+1),std::vector<size_t>(ci,ci+nnz),std::vector<double>(v,v+nnz));
#else
    sparse_matrix sparse={0};require(m_csr_create(rows,cols,nnz,rp,ci,v,&sparse));
#endif
#ifdef __cplusplus
    std::unique_ptr<linear_a::ILU0> factor;
    if(!strcmp(op,"ilu_apply")||(!strcmp(op,"gmres")&&jacobi==2))factor.reset(new linear_a::ILU0(sparse));
#else
    matrix_ilu0 factor={0};
    if(!strcmp(op,"ilu_apply")||(!strcmp(op,"gmres")&&jacobi==2))require(m_ilu0_create(&sparse,&factor));
#endif
    if(!strcmp(op,"dense")){require(m_create(rows,cols,&dense));for(size_t i=0;i<rows;i++)for(size_t p=rp[i];p<rp[i+1];p++)dense.values[i*cols+ci[p]]=v[p];}
    free(rp);free(ci);free(v);
    double checksum=0;uint64_t start=0,elapsed=0;size_t runs=iterations?iterations+3:1;
    for(size_t run=0;run<runs;run++){
        if(iterations&&run==3)start=now_ns();
        matrix out={0};
        if(!strcmp(op,"ilu_setup")){
#ifdef __cplusplus
            linear_a::ILU0 f(sparse);require(m_create(1,2,&out));out.values[0]=(double)f.size();out.values[1]=(double)f.nnz();
#else
            matrix_ilu0 f={0};require(m_ilu0_create(&sparse,&f));require(m_create(1,2,&out));out.values[0]=(double)f.factors.rows;out.values[1]=(double)f.factors.nnz;m_ilu0_free(&f);
#endif
        }else if(!strcmp(op,"ilu_apply")){
#ifdef __cplusplus
            auto x=factor->apply(std::vector<double>(b.values,b.values+count));require(m_create(rows,1,&out));if(rows)memcpy(out.values,x.data(),rows*sizeof(double));
#else
            require(m_ilu0_apply(&factor,&b,&out));
#endif
        }else if(!strcmp(op,"dense"))require(m_multiply(&dense,&b,&out));
        else if(!strcmp(op,"spmv")){
#ifdef __cplusplus
            auto result=sparse.matvec(std::vector<double>(b.values,b.values+count));require(m_create(rows,1,&out));if(rows)memcpy(out.values,result.data(),rows*sizeof(double));
#else
            require(m_csr_matvec(&sparse,&b,&out));
#endif
        }else if(!strcmp(op,"gmres")){
            size_t restart=integer(argv[12]);
#ifdef __cplusplus
            std::unique_ptr<linear_a::ILU0> current;if(jacobi==3)current.reset(new linear_a::ILU0(sparse));
            auto result=sparse.gmres(std::vector<double>(b.values,b.values+count),restart,rtol,atol,limit,jacobi==1,capture!=0,current?current.get():factor.get());
            int reason=result.reason=="converged"?0:result.reason=="iteration_limit"?1:result.reason=="breakdown"?2:result.reason=="nonfinite"?3:4;
            size_t steps=result.iterations,restart_count=result.restarts.size();const double *x=result.x.data(),*history=result.residuals.data(),*estimated=result.estimated_residuals.data();const size_t *restarts=result.restarts.data();
            std::vector<double> flat;for(const auto& f:result.iterates)flat.insert(flat.end(),f.begin(),f.end());const double *frames=flat.data();
#else
            matrix_ilu0 current={0};if(jacobi==3)require(m_ilu0_create(&sparse,&current));
            matrix_gmres_result result={0};require(m_csr_gmres_preconditioned(&sparse,b.values,count,restart,rtol,atol,limit,jacobi==1,capture!=0,jacobi==3?&current:jacobi==2?&factor:NULL,&result));m_ilu0_free(&current);
            int reason=result.reason;size_t steps=result.iterations,restart_count=result.restart_count;const double *x=result.x,*history=result.residuals,*estimated=result.estimated_residuals,*frames=result.iterates;const size_t *restarts=result.restarts;
#endif
            if(iterations){if(reason)fail("benchmark GMRES did not converge");require(m_create(rows,1,&out));if(rows)memcpy(out.values,x,rows*sizeof(double));}
            else{
                size_t length=4+rows+2*(steps+1)+restart_count+(capture?(steps+1)*rows:0);require(m_create(1,length,&out));out.values[0]=reason;out.values[1]=(double)steps;out.values[2]=(double)(steps+1);size_t offset=3;
                if(rows)memcpy(out.values+offset,x,rows*sizeof(double));offset+=rows;memcpy(out.values+offset,history,(steps+1)*sizeof(double));offset+=steps+1;memcpy(out.values+offset,estimated,(steps+1)*sizeof(double));offset+=steps+1;out.values[offset++]=(double)restart_count;
                for(size_t k=0;k<restart_count;k++)out.values[offset++]=(double)restarts[k];if(capture&&rows)memcpy(out.values+offset,frames,(steps+1)*rows*sizeof(double));
            }
#ifndef __cplusplus
            m_gmres_free(&result);
#endif
        }else{
#ifdef __cplusplus
            auto result=sparse.conjugate_gradient(std::vector<double>(b.values,b.values+count),rtol,atol,limit,jacobi!=0,capture!=0);
            int reason=result.reason=="converged"?0:result.reason=="iteration_limit"?1:result.reason=="breakdown"?2:3;
            size_t steps=result.iterations;const double *x=result.x.data(),*history=result.residuals.data();
            std::vector<double> flat;for(const auto& f:result.iterates)flat.insert(flat.end(),f.begin(),f.end());const double *frames=flat.data();
#else
            matrix_cg_result result={0};require(m_csr_cg(&sparse,b.values,count,rtol,atol,limit,jacobi!=0,capture!=0,&result));
            int reason=result.reason;size_t steps=result.iterations;const double *x=result.x,*history=result.residuals,*frames=result.iterates;
#endif
            if(iterations){if(reason)fail("benchmark CG did not converge");require(m_create(rows,1,&out));if(rows)memcpy(out.values,x,rows*sizeof(double));}
            else{
                size_t length=3+rows+steps+1+(capture?(steps+1)*rows:0);require(m_create(1,length,&out));out.values[0]=reason;out.values[1]=(double)steps;out.values[2]=(double)(steps+1);
                if(rows)memcpy(out.values+3,x,rows*sizeof(double));memcpy(out.values+3+rows,history,(steps+1)*sizeof(double));if(capture&&rows)memcpy(out.values+4+rows+steps,frames,(steps+1)*rows*sizeof(double));
            }
#ifndef __cplusplus
            m_cg_free(&result);
#endif
        }
        if(!iterations){if(strcmp(op,"cg")&&strcmp(op,"gmres")){printf("{\"rows\":1,\"cols\":%zu,\"values\":[",out.rows*out.cols);for(size_t i=0;i<out.rows*out.cols;i++)printf("%s%.17g",i?",":"",out.values[i]);puts("]}");}else{print_matrix(&out);puts("");}}
        else if(run>=3)for(size_t i=0;i<out.rows*out.cols;i++)checksum+=out.values[i];
        m_free(&out);
    }
    if(iterations){elapsed=now_ns()-start;if(!isfinite(checksum))fail("nonfinite checksum");printf("{\"elapsed_ns\":%" PRIu64 ",\"iterations\":%zu,\"checksum\":%.17g}\n",elapsed,iterations,checksum);}
#ifndef __cplusplus
    m_ilu0_free(&factor);m_csr_free(&sparse);
#endif
    m_free(&b);m_free(&dense);return 0;
}
