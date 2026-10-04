#ifdef __cplusplus
#include <memory>
#endif
/* Shared protocol driver; native C++ calls its public CSRMatrix interface. */
static double sparse_token(void){double x;if(scanf("%lf",&x)!=1||!isfinite(x))fail("invalid sparse input token");return x;}
static size_t sparse_index(void){double x=sparse_token();if(x<0||x>=(double)SIZE_MAX||floor(x)!=x)fail("invalid CSR index");return (size_t)x;}
static int sparse_cli(int argc,char **argv){
    if(argc!=(argc>2&&!strcmp(argv[2],"gmres")?13:12))fail("expected sparse OP ROWS COLS NNZ ITERATIONS RTOL ATOL LIMIT JACOBI CAPTURE");
    const char *op=argv[2];if(strcmp(op,"spmv")&&strcmp(op,"dense")&&strcmp(op,"cg")&&strcmp(op,"gmres")&&strcmp(op,"ic0_factor")&&strcmp(op,"ic0_apply")&&strcmp(op,"ilu_setup")&&strcmp(op,"ilu_apply")&&strcmp(op,"rcm")&&strcmp(op,"amd")&&strcmp(op,"permute")&&strcmp(op,"permutation_check")&&strcmp(op,"rcm_solve")&&strcmp(op,"ilu_solve")&&strcmp(op,"chol_symbolic")&&strcmp(op,"chol_factor")&&strcmp(op,"chol_solve")&&strcmp(op,"chol_total")&&strcmp(op,"chol_rcm_total")&&strcmp(op,"chol_amd_total"))fail("invalid sparse operation");
    size_t rows=integer(argv[3]),cols=integer(argv[4]),nnz=integer(argv[5]),iterations=integer(argv[6]),limit=integer(argv[9]),jacobi=integer(argv[10]),capture=integer(argv[11]);
    double rtol=number(argv[7]),atol=number(argv[8]);if(jacobi>3||(strcmp(op,"gmres")&&strcmp(op,"cg")&&jacobi>1)||capture>1||rows==SIZE_MAX||rows>SIZE_MAX/sizeof(size_t)-1||nnz>SIZE_MAX/sizeof(double))fail("invalid sparse dimensions/options");
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
    std::unique_ptr<linear_a::IC0> ic;if(!strcmp(op,"ic0_apply")||(!strcmp(op,"cg")&&jacobi==2))ic.reset(new linear_a::IC0(sparse));
    std::unique_ptr<linear_a::ILU0> factor;
    if(!strcmp(op,"ilu_apply")||(!strcmp(op,"gmres")&&jacobi==2))factor.reset(new linear_a::ILU0(sparse));
#else
    matrix_ic0 ic={0};if(!strcmp(op,"ic0_apply")||(!strcmp(op,"cg")&&jacobi==2))require(m_ic0_create(&sparse,&ic));
    matrix_ilu0 factor={0};
    if(!strcmp(op,"ilu_apply")||(!strcmp(op,"gmres")&&jacobi==2))require(m_ilu0_create(&sparse,&factor));
#endif
#ifdef __cplusplus
    std::unique_ptr<linear_a::SparseCholeskySymbolic> plan;
    std::unique_ptr<linear_a::SparseCholesky> chol;
    if(!strcmp(op,"chol_factor")||!strcmp(op,"chol_solve"))plan.reset(new linear_a::SparseCholeskySymbolic(sparse));
    if(!strcmp(op,"chol_solve"))chol.reset(new linear_a::SparseCholesky(plan->factorize(sparse)));
#else
    matrix_cholesky_symbolic plan={0};matrix_cholesky chol={0};
    if(!strcmp(op,"chol_factor")||!strcmp(op,"chol_solve"))require(m_cholesky_analyze(&sparse,&plan));
    if(!strcmp(op,"chol_solve"))require(m_cholesky_factorize(&plan,&sparse,&chol));
#endif
    if(!strcmp(op,"dense")){require(m_create(rows,cols,&dense));for(size_t i=0;i<rows;i++)for(size_t p=rp[i];p<rp[i+1];p++)dense.values[i*cols+ci[p]]=v[p];}
    free(rp);free(ci);free(v);
    double checksum=0;uint64_t start=0,elapsed=0;size_t runs=iterations?iterations+3:1;
    for(size_t run=0;run<runs;run++){
        if(iterations&&run==3)start=now_ns();
        matrix out={0};
        if(!strcmp(op,"ic0_factor")){
#ifdef __cplusplus
            linear_a::IC0 f(sparse);auto l=f.lower();size_t m=l.nnz();require(m_create(1,rows+1+2*m,&out));size_t k=0;for(auto v:l.row_offsets())out.values[k++]=(double)v;for(auto v:l.column_indices())out.values[k++]=(double)v;for(auto v:l.values())out.values[k++]=v;
#else
            matrix_ic0 f={0};require(m_ic0_create(&sparse,&f));size_t m=f.lower.nnz;require(m_create(1,rows+1+2*m,&out));size_t k=0;for(size_t i=0;i<=rows;i++)out.values[k++]=(double)f.lower.offsets[i];for(size_t p=0;p<m;p++)out.values[k++]=(double)f.lower.indices[p];for(size_t p=0;p<m;p++)out.values[k++]=f.lower.values[p];m_ic0_free(&f);
#endif
        }else if(!strcmp(op,"ic0_apply")){
#ifdef __cplusplus
            auto x=ic->apply(std::vector<double>(b.values,b.values+count));require(m_create(1,x.size(),&out));if(x.size())memcpy(out.values,x.data(),x.size()*sizeof(double));
#else
            require(m_ic0_apply(&ic,&b,&out));
#endif
        }else if(!strcmp(op,"chol_symbolic")){
#ifdef __cplusplus
            linear_a::SparseCholeskySymbolic s(sparse);auto rp=s.row_offsets(),ci=s.column_indices();auto steps=s.fill_steps();size_t m=s.nnz();require(m_create(1,rows+1+2*m,&out));size_t k=0;for(auto x:rp)out.values[k++]=(double)x;for(auto x:ci)out.values[k++]=(double)x;for(auto x:steps)out.values[k++]=(double)x;
#else
            matrix_cholesky_symbolic s={0};require(m_cholesky_analyze(&sparse,&s));size_t m=s.pattern.nnz;require(m_create(1,rows+1+2*m,&out));size_t k=0;for(size_t i=0;i<=rows;i++)out.values[k++]=(double)s.pattern.offsets[i];for(size_t p=0;p<m;p++)out.values[k++]=(double)s.pattern.indices[p];for(size_t p=0;p<m;p++)out.values[k++]=s.steps[p]==SIZE_MAX?-1:(double)s.steps[p];m_cholesky_symbolic_free(&s);
#endif
        }else if(!strcmp(op,"chol_factor")){
#ifdef __cplusplus
            auto f=plan->factorize(sparse);auto l=f.lower();size_t m=l.nnz();require(m_create(1,rows+1+2*m,&out));size_t k=0;for(auto x:l.row_offsets())out.values[k++]=(double)x;for(auto x:l.column_indices())out.values[k++]=(double)x;for(auto x:l.values())out.values[k++]=x;
#else
            matrix_cholesky f={0};require(m_cholesky_factorize(&plan,&sparse,&f));size_t m=f.lower.nnz;require(m_create(1,rows+1+2*m,&out));size_t k=0;for(size_t i=0;i<=rows;i++)out.values[k++]=(double)f.lower.offsets[i];for(size_t p=0;p<m;p++)out.values[k++]=(double)f.lower.indices[p];for(size_t p=0;p<m;p++)out.values[k++]=f.lower.values[p];m_cholesky_free(&f);
#endif
        }else if(!strcmp(op,"chol_solve")){
#ifdef __cplusplus
            auto x=chol->solve(std::vector<double>(b.values,b.values+count));require(m_create(rows,1,&out));for(size_t i=0;i<rows;i++)out.values[i]=x[i];
#else
            require(m_cholesky_solve(&chol,&b,&out));
#endif
        }else if(!strcmp(op,"chol_total")||!strcmp(op,"chol_rcm_total")||!strcmp(op,"chol_amd_total")){
            int amd=!strcmp(op,"chol_amd_total"),reorder=amd||!strcmp(op,"chol_rcm_total");
#ifdef __cplusplus
            std::unique_ptr<linear_a::CSRMatrix> reordered;std::vector<size_t> p;std::vector<double> rhs(b.values,b.values+count);
            if(reorder){p=amd?sparse.approximate_minimum_degree():sparse.reverse_cuthill_mckee();reordered.reset(new linear_a::CSRMatrix(sparse.permute_symmetric(p)));rhs=linear_a::CSRMatrix::permute_vector(p,rhs);}
            const auto& q=reordered?*reordered:sparse;linear_a::SparseCholeskySymbolic s(q);auto f=s.factorize(q);auto x=f.solve(rhs);if(reorder)x=linear_a::CSRMatrix::permute_vector(p,x,true);require(m_create(rows,1,&out));for(size_t i=0;i<rows;i++)out.values[i]=x[i];
#else
            sparse_matrix reordered={0};const sparse_matrix *q=&sparse;size_t *p=NULL;matrix rhs={0};const matrix *right=&b;
            if(reorder){p=(size_t*)calloc(rows?rows:1,sizeof(size_t));if(!p)fail("allocation failed");require(amd?m_csr_amd(&sparse,p,rows):m_csr_rcm(&sparse,p,rows));require(m_csr_permute(&sparse,p,rows,&reordered));q=&reordered;require(m_create(rows,1,&rhs));require(m_permute_vector(p,rows,b.values,false,rhs.values));right=&rhs;}
            matrix_cholesky_symbolic s={0};matrix_cholesky f={0};matrix x={0};require(m_cholesky_analyze(q,&s));require(m_cholesky_factorize(&s,q,&f));require(m_cholesky_solve(&f,right,&x));require(m_create(rows,1,&out));if(reorder)require(m_permute_vector(p,rows,x.values,true,out.values));else if(rows)memcpy(out.values,x.values,rows*sizeof(double));m_free(&x);m_cholesky_free(&f);m_cholesky_symbolic_free(&s);m_csr_free(&reordered);m_free(&rhs);free(p);
#endif
        }else if(!strcmp(op,"rcm_solve")||!strcmp(op,"ilu_solve")){
            int reorder=!strcmp(op,"rcm_solve");require(m_create(1,rows+1,&out));
#ifdef __cplusplus
            std::unique_ptr<linear_a::CSRMatrix> reordered;std::vector<size_t> p;std::vector<double> rhs(b.values,b.values+count);
            if(reorder){p=sparse.reverse_cuthill_mckee();reordered.reset(new linear_a::CSRMatrix(sparse.permute_symmetric(p)));rhs=linear_a::CSRMatrix::permute_vector(p,rhs);}
            const auto& q=reordered?*reordered:sparse;linear_a::ILU0 f(q);auto result=q.gmres(rhs,20,rtol,atol,limit,false,false,&f);if(!result.converged)fail("ordering solve failed");auto x=reorder?linear_a::CSRMatrix::permute_vector(p,result.x,true):result.x;out.values[0]=(double)result.iterations;for(size_t i=0;i<rows;i++)out.values[i+1]=x[i];
#else
            sparse_matrix reordered={0};const sparse_matrix *q=&sparse;size_t *p=NULL;matrix rhs={0};const double *right=b.values;
            if(reorder){p=(size_t*)calloc(rows?rows:1,sizeof(size_t));if(!p)fail("allocation failed");require(m_csr_rcm(&sparse,p,rows));require(m_csr_permute(&sparse,p,rows,&reordered));q=&reordered;require(m_create(rows,1,&rhs));require(m_permute_vector(p,rows,b.values,false,rhs.values));right=rhs.values;}
            matrix_ilu0 f={0};require(m_ilu0_create(q,&f));matrix_gmres_result result={0};require(m_csr_gmres_preconditioned(q,right,count,20,rtol,atol,limit,false,false,&f,&result));if(result.reason)fail("ordering solve failed");out.values[0]=(double)result.iterations;if(reorder)require(m_permute_vector(p,rows,result.x,true,out.values+1));else if(rows)memcpy(out.values+1,result.x,rows*sizeof(double));m_gmres_free(&result);m_ilu0_free(&f);m_csr_free(&reordered);m_free(&rhs);free(p);
#endif
        }else if(!strcmp(op,"amd")){
            require(m_create(rows,1,&out));
#ifdef __cplusplus
            auto p=sparse.approximate_minimum_degree();for(size_t i=0;i<rows;i++)out.values[i]=(double)p[i];
#else
            size_t *p=(size_t*)calloc(rows?rows:1,sizeof(size_t));if(!p)fail("allocation failed");require(m_csr_amd(&sparse,p,rows));for(size_t i=0;i<rows;i++)out.values[i]=(double)p[i];free(p);
#endif
        }else if(!strcmp(op,"rcm")){
            require(m_create(rows,1,&out));
#ifdef __cplusplus
            auto p=sparse.reverse_cuthill_mckee();for(size_t i=0;i<rows;i++)out.values[i]=(double)p[i];
#else
            size_t *p=(size_t*)calloc(rows?rows:1,sizeof(size_t));if(!p)fail("allocation failed");require(m_csr_rcm(&sparse,p,rows));for(size_t i=0;i<rows;i++)out.values[i]=(double)p[i];free(p);
#endif
        }else if(!strcmp(op,"permute")||!strcmp(op,"permutation_check")){
            size_t *p=(size_t*)calloc(count?count:1,sizeof(size_t));if(!p)fail("allocation failed");for(size_t i=0;i<count;i++){double x=b.values[i];if(x<0||x>=(double)rows||floor(x)!=x)fail("invalid permutation");p[i]=(size_t)x;}
#ifdef __cplusplus
            std::vector<size_t> order(p,p+count);auto q=sparse.permute_symmetric(order);const auto& qr=q.row_offsets();const auto& qc=q.column_indices();const auto& qv=q.values();
#else
            sparse_matrix q={0};require(m_csr_permute(&sparse,p,count,&q));const size_t *qr=q.offsets,*qc=q.indices;const double *qv=q.values;
#endif
            if(!strcmp(op,"permute")){require(m_create(1,rows+1+2*nnz,&out));size_t k=0;for(size_t i=0;i<=rows;i++)out.values[k++]=(double)qr[i];for(size_t i=0;i<nnz;i++)out.values[k++]=(double)qc[i];for(size_t i=0;i<nnz;i++)out.values[k++]=qv[i];}
            else{
                require(m_create(1,3*rows,&out));
#ifdef __cplusplus
                std::vector<double> x(rows);for(size_t i=0;i<rows;i++)x[i]=(double)(i+1);auto y=linear_a::CSRMatrix::permute_vector(order,x);auto back=linear_a::CSRMatrix::permute_vector(order,y,true);auto z=linear_a::CSRMatrix::permute_vector(order,q.matvec(y),true);for(size_t i=0;i<rows;i++){out.values[i]=y[i];out.values[rows+i]=back[i];out.values[2*rows+i]=z[i];}
#else
                matrix x={0},y={0},z={0};require(m_create(rows,1,&x));require(m_create(rows,1,&y));for(size_t i=0;i<rows;i++)x.values[i]=(double)(i+1);require(m_permute_vector(p,rows,x.values,false,y.values));if(rows)memcpy(out.values,y.values,rows*sizeof(double));require(m_permute_vector(p,rows,y.values,true,out.values+rows));require(m_csr_matvec(&q,&y,&z));require(m_permute_vector(p,rows,z.values,true,out.values+2*rows));m_free(&x);m_free(&y);m_free(&z);
#endif
            }
#ifndef __cplusplus
            m_csr_free(&q);
#endif
            free(p);
        }else if(!strcmp(op,"ilu_setup")){
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
            std::unique_ptr<linear_a::IC0> current;if(jacobi==3)current.reset(new linear_a::IC0(sparse));
            auto result=sparse.conjugate_gradient(std::vector<double>(b.values,b.values+count),rtol,atol,limit,jacobi==1,capture!=0,jacobi==3?current.get():ic.get());
            int reason=result.reason=="converged"?0:result.reason=="iteration_limit"?1:result.reason=="breakdown"?2:3;
            size_t steps=result.iterations;const double *x=result.x.data(),*history=result.residuals.data();
            std::vector<double> flat;for(const auto& f:result.iterates)flat.insert(flat.end(),f.begin(),f.end());const double *frames=flat.data();
#else
            matrix_ic0 current={0};if(jacobi==3)require(m_ic0_create(&sparse,&current));
            matrix_cg_result result={0};require(m_csr_cg_preconditioned(&sparse,b.values,count,rtol,atol,limit,jacobi==1,capture!=0,jacobi==3?&current:jacobi==2?&ic:NULL,&result));m_ic0_free(&current);
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
    m_ic0_free(&ic);m_cholesky_free(&chol);m_cholesky_symbolic_free(&plan);m_ilu0_free(&factor);m_csr_free(&sparse);
#endif
    m_free(&b);m_free(&dense);return 0;
}
