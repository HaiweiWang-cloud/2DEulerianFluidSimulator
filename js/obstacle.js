class Obstacle {
    constructor(x, y, r, pxsPerLength) {
        this.x = x;
        this.y = y;
        this.r = r;
        this.vx = 0;
        this.vy = 0;
        this.pxsPerLength = pxsPerLength
    }

    update(dt) {
        this.x += this.vx * dt;
        this.y += this.vy * dt;
    }

    setObstacle(fluid) {
        const n = fluid.Ny;
    
        for (let i = 1; i < fluid.Nx-2; i++) {
            for (let j = 1; j < fluid.Ny-2; j++) {
    
                const dx = (i + 0.5) * fluid.h - this.x;
                const dy = (j + 0.5) * fluid.h - this.y;
                
                if (dx * dx + dy * dy < this.r * this.r) {
                    fluid.s[i*n + j] = 0.0;
                    fluid.u[i*n + j] = this.vx;
                    fluid.u[(i+1)*n + j] = this.vx;
                    fluid.v[i*n + j] = this.vy;
                    fluid.v[i*n + j] = this.vy;
                } else {
                    fluid.s[i*n+j] = 1.0;
                }
                
            }
        }
    }

    draw(ctx) {
        ctx.fillStyle = "black";
        ctx.beginPath();
        ctx.arc(this.x * this.pxsPerLength, this.y * this.pxsPerLength, this.r * this.pxsPerLength, 0, 2 * Math.PI);
        ctx.fill();
    }
}

