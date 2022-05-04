namespace Harmonic.methods;
public abstract class LUMethods
{   
    private SLAU slau; 
    protected Vector lgl, lgu, ldi;
    protected int N, N_el;

    //* Инициализация 
    protected void Init(SLAU slau) {
        this.slau = slau;
        this.N    = slau.N;
        this.N_el = slau.N_el;
        lgl       = new Vector(N_el);
        lgu       = new Vector(N_el);
        ldi       = new Vector(N);
        Vector.Copy(slau.gl, lgl);
        Vector.Copy(slau.gu, lgu);
        Vector.Copy(slau.di, ldi);
    }

    //* Инициализация Profile
    protected void Init_profile(SLAU slau) {
        this.N    = slau.N;
    
        int[] newig = new int[this.N + 1];
        for (int i = 1; i <= this.N; i++) {
            if ((slau.ig[i] - slau.ig[i - 1]) > 0)
                newig[i] = newig[i - 1] + (i - slau.jg[slau.ig[i - 1]]);
            else 
                newig[i] = newig[i - 1];
        }

        slau.N_el = newig[this.N];
        this.N_el = slau.N_el;

        Vector newgl = new Vector(this.N_el);
        Vector newgu = new Vector(this.N_el);

        for (int i = 0; i < this.N; i++) {
            int col = i - (newig[i + 1] - newig[i]);
            int p = slau.ig[i];
            for(int j = newig[i]; j < newig[i + 1]; j++, col++){
                if(col == slau.jg[p]){
                    newgu[j] = slau.gu[p];
                    newgl[j] = slau.gl[p];
                    p++;
                }
                else
                    newgu[j] = newgl[j] = 0;
            }
        }

        slau.ig = new int[this.N + 1];
        slau.gl = new Vector(this.N_el);
        slau.gu = new Vector(this.N_el);
        Array.Copy(newig, slau.ig, newig.Length);
        Vector.Copy(newgl, slau.gl);
        Vector.Copy(newgu, slau.gu);
        
        this.slau = slau;
    }

    //* lU-разложение 
    protected void LU_decomposition() {
        double sumDI, sumGL, sumGU;
        for (int i = 0; i < N; i++) {
            sumDI = 0;
            for (int j = slau.ig[i]; j < slau.ig[i + 1]; j++) {
                sumGL = sumGU = 0;
                int p_s = slau.ig[slau.jg[j]], p_s1 = slau.ig[slau.jg[j] + 1];
                for (int k = slau.ig[i]; k < j; k++)
                    for (int p = p_s; p < p_s1; p++)
                        if (slau.jg[k] == slau.jg[p]) {
                            sumGL += lgl[k]*lgu[p];
                            sumGU += lgl[p]*lgu[k];
                            p_s++;
                        }
                lgl[j] -= sumGL;
                lgu[j]  = (lgu[j] - sumGU) / ldi[slau.jg[j]]; 
                sumDI  += lgl[j]*lgu[j];
            }
            ldi[i] -= sumDI;
        }
    }

    //* Решение СЛАУ L^(-1)
    protected Vector Solve_L(Vector f) {
        Vector x = new Vector(N);
        double sum;
        for (int i = 0; i < N; i++) {
            sum = 0;
            for (int j = slau.ig[i]; j < slau.ig[i + 1]; j++)
                sum += lgl[j]*x[slau.jg[j]];
            x[i] = (f[i] - sum) / ldi[i];
        }
        return x;
    }

    //* Решение СЛАУ U^(-1)
    protected Vector Solve_U(Vector f) {
        Vector x = new Vector(N);
        Vector f_copy = new Vector(N);
        Vector.Copy(f, f_copy);
        for (int i = N - 1; i >= 0; i--) {
            x[i] = f_copy[i] / ldi[i];
            for (int j = slau.ig[i]; j < slau.ig[i + 1]; j++)
                f_copy[slau.jg[j]] -= lgu[j]*x[i];
        }
        return x;
    }


    //* lU-разложение Profile
    protected void LU_decomposition_profile() {
        for(int i = 0; i < slau.N; i++){
            int i0 = slau.ig[i]; 
            int i1 = slau.ig[i + 1];
            int j = i - (i1 - i0);
            double sd = 0;
            for(int m = i0; m < i1; m++,j++){
                double sl = 0;
                double su = 0;
                int j0 = j < 0 ? -33686019 : slau.ig[j]; 
                int j1 = slau.ig[j + 1];
                int mi = i0;
                int mj = j0;
                int kol_i = m - i0;
                int kol_j = j1 - j0;
                int kol_r = kol_i - kol_j;
                if(kol_r < 0) mj -= kol_r;
                else mi += kol_r;
                for(; mi<m; mi++, mj++){
                    sl += slau.gl[mi]*slau.gu[mj];
                    su += slau.gu[mi]*slau.gl[mj];
                }
                slau.gl[m] = slau.gl[m] - sl;
                slau.gu[m] = j < 0 ? 0 : (slau.gu[m] - su) / slau.di[j];
                sd += slau.gl[m]*slau.gu[m];
            }
            slau.di[i] = slau.di[i] - sd;
        }
    }

    //* Прямой ход
    protected Vector normalRunning(Vector f) {     
        Vector y = new Vector(slau.N);
        for (int i = 0; i < slau.N; i++) {
            int i0 = slau.ig[i], 
                i1 = slau.ig[i + 1], 
                j = i - (i1 - i0);
            double sum = 0;
            for (int k = i0; k < i1; k++, j++)
                sum += j < 0 ? 0 : y[j] * slau.gl[k];
            y[i] = (f[i] - sum) / slau.di[i];
        }
        return y;
    }

    //* Обратный ход
    protected Vector reverseRunning(Vector y) {      
        Vector x = new Vector(slau.N);
        for (int i = slau.N - 1; i >= 0; i--) {
            int i0 = slau.ig[i], 
                i1 = slau.ig[i + 1],
                j = i - (i1 - i0);
            double xi = y[i];
            for (int k = i0; k < i1; k++, j++) {
                if (j >= 0)
                    y[j] -= xi * slau.gu[k];
            }      
            x[i] = xi;
        }
        return x;
    }
}
