function jet(max, min, value) {
    let val = Math.min(Math.max(value, min), max-0.0001);
    const d = max - min;
    val = d == 0 ? 0.5 : (val - min) / d;
    const m = 0.25;
    const num = Math.floor(val / m);
    const s = (val - num * m) / m;
    let r, g, b;

    switch (num) {
        case 0: r = 0; g=s; b = 1; break;
        case 1: r = 0; g=1; b = 1-s; break;
        case 2: r=s; g=1; b=0; break;
        case 3: r=1; g=1-s; b=0; break;
    }

    return [r*255, g*255, b*255, 255];
}

function grayscale(max, min, value) {
    let val = Math.min(Math.max(value, min), max-0.0001);
    const d = max - min;
    val = d == 0 ? 0.5 : (val - min) / d;
    
    return [255*val, 255*val, 255*val, 255];
}