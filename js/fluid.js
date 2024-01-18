const U_FIELD = 1;
const V_FIELD = 2;
const S_FIELD = 0;

class Fluid {
    constructor(canvas, density, N, M, h, gridScale) {
        this.canvas = canvas;
        this.ctx = this.canvas.getContext("2d");
        this.gridScale = Math.floor(gridScale);
        this.rho = density;
        this.Nx = N;
        this.Ny = M;
        this.h = h;
        this.overRelaxation = 1.9;
        this.numCells = this.Nx * this.Ny;
        this.u = new Float32Array(this.numCells); // x velocity
        this.v = new Float32Array(this.numCells); // y velocity
        this.p = new Float32Array(this.numCells); // pressure
        this.s = new Float32Array(this.numCells) // boundaries array
        this.m = new Float32Array(this.numCells) // smoke array
        this.newU = new Float32Array(this.numCells); // helpers for updating
        this.newV = new Float32Array(this.numCells);
        this.newM = new Float32Array(this.numCells);
        this.m.fill(1.0);
    }

    // updating dynamics
    update(dt, gravity, numIters) {
        this.#integrate(dt, gravity);
        this.p.fill(0.0);
        this.#solveIncompressibility(numIters, dt);

        this.#extrapolate();
        this.#advectVelocity(dt);
        this.#advectSmoke(dt);
    }

    #integrate(dt, gravity) {
        const n = this.Ny
        for (let i=1; i < this.Nx; i++) {
            for (let j=1; j < this.Ny-1; j++) {
                if (this.s[i * n + j] != 0.0 && this.s[i*n + j - 1] != 0.0) {
                    this.v[i * n+ j] += gravity * dt;
                }
            }
        }
    }

    #solveIncompressibility(numIters, dt) {
        const n = this.Ny;
        const cp = this.rho * this.h / dt;
        const o = this.overRelaxation;

        for (let iter = 0; iter < numIters; iter++) {

            for (let i=1; i<this.Nx-1; i++) {
                for (let j=1; j<this.Ny-1; j++) {
                    
                    if (this.s[i * n + j] == 0.0) {
                        continue;
                    }
                        
                    let s = this.s[i*n + j];
                    const sx0 = this.s[(i-1) * n + j];
                    const sx1 = this.s[(i+1) * n + j];
                    const sy0 = this.s[i * n + j-1];
                    const sy1 = this.s[i * n + j+1];
                    s = sx0 + sx1 + sy0 + sy1;

                    if (s==0.0) { 
                        continue;
                    }

                    const div = this.u[(i+1)*n + j] - this.u[i*n + j] + this.v[i*n + j+1] - this.v[i*n + j];
                   
                    let p = -div / s * o;
                    this.p[i*n + j] += cp * p;
                    
                    this.u[i*n+j] -= sx0 * p;
                    this.u[(i+1)*n+j] += sx1 * p;
                    this.v[i*n+j] -= sy0 * p;
                    this.v[i*n+j+1] += sy1 * p;
                }
            }
        }
    }

    #extrapolate() {
        // Matched boundary conditions
        const n = this.Ny;
        for (let i = 0; i < this.Nx; i++) {
            this.u[i*n + 0] = this.u[i*n + 1];
            this.u[i*n+this.Ny-1] = this.u[i*n+this.Ny-2];
        }
        for (let j=0; j < this.Ny; j++) {
            this.v[0*n + j] = this.v[1*n+j];
            this.v[(this.Nx-1)*n + j] = this.v[(this.Nx-2)*n + j];
        }
    }

    #avgV(i, j) {
        const n = this.Ny;
        return 0.25 * (this.v[i*n + j] + this.v[(i-1)*n+j] + this.v[(i-1)*n+j+1] + this.v[i*n+j+1])
    }

    #avgU(i, j) {
        const n = this.Ny;
        return 0.25 * (this.u[i*n + j] + this.u[(i+1)*n+j] + this.u[(i+1)*n+j-1] + this.u[i*n+j-1])
    }

    sampleField(xc, yc, field) {
        const n = this.Ny;
        const h = this.h;
        const h1 = 1.0 / h;
        const h2 = 0.5 * h;

        const x = Math.max(Math.min(xc, this.Nx * h), h);
        const y = Math.max(Math.min(yc, this.Ny * h), h);

        let dx = 0.0;
        let dy = 0.0;

        let f;

        switch (field) {
            case U_FIELD: f = this.u; dy = h2; break;
            case V_FIELD: f = this.v; dx = h2; break;
            case S_FIELD: f = this.m; dx = h2; dy = h2; break;
        }

        const x0 = Math.min(Math.floor((x-dx)*h1), this.Nx-1);
        const tx = ((x-dx) - x0*h) * h1;
        const x1 = Math.min(x0 + 1, this.Nx-1);

        const y0 = Math.min(Math.floor((y-dy)*h1), this.Ny-1);
        const ty = ((y-dy) - y0*h) * h1;
        const y1 = Math.min(y0 + 1, this.Ny-1);

        const sx = 1.0 - tx;
        const sy = 1.0 - ty;

        return sx * sy * f[x0*n + y0] + tx*sy * f[x1*n + y0] + tx*ty * f[x1*n+y1] + sx*ty * f[x0*n+y1];
    }

    #advectVelocity(dt) {
        this.newU.set(this.u);
        this.newV.set(this.v);

        const n = this.Ny;
        const h = this.h;
        const h2 = 0.5 * h;

        for (let i = 1; i < this.Nx; i++) {
            for (let j = 1; j < this.Ny; j++) {

                // u component
                if (this.s[i*n+j] != 0.0 && this.s[(i-1)*n+j] != 0.0 && j < this.Ny-1) {
                    // find the fluid particle at cell i, j
                    let x = i * h;
                    let y = j * h + h2;
                    let u = this.u[i*n+j];
                    let v = this.#avgV(i, j);
                    // find particle which moves to the cell
                    x = x - dt * u;
                    y = y - dt * v;
                    u = this.sampleField(x, y, U_FIELD);
                    this.newU[i*n + j] = u;
                }
                
                // v component
                if (this.s[i*n+j] != 0.0 && this.s[i*n+j-1] != 0 && i < this.Nx - 1) {
                    let x = i * h + h2;
                    let y = j * h;
                    let u = this.#avgU(i, j);
                    let v = this.v[i*n + j];
                    x = x - dt * u;
                    y = y - dt * v;
                    v = this.sampleField(x, y, V_FIELD);
                    this.newV[i*n + j] = v;
                }
            }
        }
        this.u.set(this.newU);
        this.v.set(this.newV);
    }

    #advectSmoke(dt) {
        this.newM.set(this.m);

        const n = this.Ny;
        const h = this.h;
        const h2 = 0.5 * h;

        for (let i = 1; i < this.Nx-1; i++) {
            for (let j = 1; j < this.Ny-1; j++) {

                if (this.s[i*n + j] != 0.0) {
                    const u = (this.u[i*n+j] + this.u[(i+1)*n + j]) * 0.5;
                    const v = (this.v[i*n +j] + this.v[i*n + j+1]) * 0.5;
                    const x = i*h + h2 - dt*u;
                    const y = j*h + h2 - dt*v;

                    this.newM[i*n + j] = this.sampleField(x,y, S_FIELD);
                }
            }
        }
        this.m.set(this.newM);
    }

    // draw methods

    draw({pressure = true, smoke = false, streamline = false, velocity = false} = {}) {
        
        // 0: pressure
        if (pressure) {
            this.#drawPressure();
        } else if (smoke) {
            this.#drawSmoke();
        }

        if (streamline) {
            this.#drawStreamline();
        } else if (velocity) {
            this.#drawVelocity();
        }
    }

    #drawScalarFieldColor(max, min, scalarField, colormap) {
        const id = this.ctx.getImageData(0, 0, this.canvas.width, this.canvas.height);

        for (let i=0; i<this.Nx; i++) {
            for (let j=0; j<this.Ny; j++) {
                const color = colormap(max, min, scalarField[i*this.Ny + j]);
                
                const r = color[0];
                const g = color[1];
                const b = color[2];
                const a = color[3];

                const x = Math.floor(i * this.gridScale);
                const y = Math.floor(j * this.gridScale);
                for (let yi = y; yi < y + this.gridScale; yi++) {
                    let p = 4 * (yi * this.canvas.width + x);
                    for (let xi = x; xi < x + this.gridScale; xi++){
                        id.data[p++] = r;
                        id.data[p++] = g;
                        id.data[p++] = b;
                        id.data[p++] = a;
                    }
                }
            }
        }

        this.ctx.putImageData(id, 0, 0);
    }

    #drawPressure() {
        // Find max and min values
        let minP = this.p[0];
        let maxP = this.p[0];
        for (let i=0; i < this.numCells; i++) {
            minP = Math.min(minP, this.p[i]);
            maxP = Math.max(maxP, this.p[i]);
        }

        this.#drawScalarFieldColor(maxP, minP, this.p, jet);

    }

    #drawSmoke() {
        let minM = this.m[0];
        let maxM = this.m[0];
        for (let i=0; i < this.numCells; i++) {
            minM = Math.min(minM, this.m[i]);
            maxM = Math.max(maxM, this.m[i]);
        }
        this.#drawScalarFieldColor(maxM, minM, this.m, grayscale);
    }

    #drawStreamline() {
        return
    }

    #drawVelocity() {
        return
    }
}