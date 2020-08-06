

/**
 * Clone an object
 * @param {*} obj 
 * @returns {*} a copy of the object
 */
function clone(obj) {
    if (typeof obj !== 'object') return obj;
    else if (obj === null) return null;

    let a = (Array.isArray(obj)) ? [] : {};
    for (key in obj) a[key] = clone(obj[key]);

    return a;
}


/**
 * Represents a vector, used for vector maths
 */
export class Vector {
    v = []
    order = 0;

    /**
     * Make a vector
     * @param  {number[]} v 
     */
    constructor(v) {
        this.v = v;
        this.order = v.length;
    }

    toArray() {
        let arr = [];
        this.v.forEach((n) => { arr.push(n); });
        return arr;
    }

    /**
     * Convert to a single column matirx
     * @returns {Matrix}
     */
    toMatrix() {
        let m = [];
        this.v.forEach((n) => { m.push([n]); });
        return new Matrix(m);
    }

    /**
     * Add two vectors together
     * @param {Vector} w
     * @returns {vector} 
     */
    add(w) {
        let arr = [];
        if (Array.isArray(w)) arr = w;
        else if (w instanceof Vector) arr = w.toArray();
        else throw 'w must be a vector or array';

        if (arr.length != this.v.length) throw 'w must be of the same order as this';

        let res = [];
        for (let i = 0; i < arr.length; i++) res.push(this.v[i] + arr[i]);

        return new Vector(res);
    }

    /**
     * Subtract w from this
     * @param {Vector} w
     * @returns {vector} 
     */
    sub(w) {
        let arr = [];
        if (Array.isArray(w)) arr = w;
        else if (w instanceof Vector) arr = w.toArray();
        else throw 'w must be a vector or array';

        if (arr.length != this.v.length) throw 'w must be of the same order as this';

        let res = [];
        for (let i = 0; i < arr.length; i++) res.push(this.v[i] - arr[i]);

        return new Vector(res);
    }

    /**
     * Compute the scalar product of this and w
     * @param {Vector} w 
     * @return {number}
     */
    scalarProd(w) {
        let arr = [];
        if (Array.isArray(w)) arr = w;
        else if (w instanceof Vector) arr = w.toArray();
        else throw 'w must be a vector or array';

        if (arr.length != this.v.length) throw 'w must be of the same order as this';

        let res = 0;
        for (let i = 0; i < arr.length; i++) res += (this.v[i] * arr[i]);

        return res;
    }

    /**
     * Compute the vector or cross product of this and w
     * @param {Vector} w a 3 vector
     * @return {Vector}
     */
    vectorProd(w) {
        let arr = [];
        if (Array.isArray(w)) arr = w;
        else if (w instanceof Vector) arr = w.toArray();
        else throw 'w must be a vector or array';

        if (arr.length != 3 || this.v.length != 3) throw 'w and this must be of order 3';

        return new Vector([
            this.v[1]*arr[2] - this.v[2]*arr[1],
            this.v[0]*arr[2] - this.v[2]*arr[0],
            this.v[0]*arr[1] - this.v[1]*arr[0]
        ]);
    }

    /**
     * Compute the inner product
     * for euclidian vectors is equivalent to the 
     * dot (or scalar) prodcut, hence this is just
     * a wrapper
     * @param {Vector} w 
     * @returns {number}
     */
    innerProd(w) {
        return this.scalarProd(w);
    }

    /**
     * Calculate the magnitude of the vector
     * @returns {number}
     */
    abs() {
        let sum = 0;
        this.v.forEach((n) => { sum += Math.pow(n, 2); })
        return Math.sqrt(sum);
    }

    /**
     * Multiply this by a constant
     * @param {number} c 
     * @returns {Vector}
     */
    mul(c) {
        let res = [];
        this.v.forEach((n) => { res.push(n*c); });
        return new Vector(res);
    }

    /**
     * Divide this by a constant
     * @param {number} c 
     * @returns {Vector}
     */
    div(c) {
        let res = [];
        this.v.forEach((n) => { res.push(n/c); });
        return new Vector(res);
    }

    /**
     * Check if the vector is a zero vector
     * @returns {boolean}
     */
    isZero() {
        let t = true;
        this.v.forEach((n) => { if (t != 0) t = false; });
        return t;
    }

    /**
     * clone a vector
     * @param {Vector} v 
     */
    static clone(v) {
        return new Vector(v.toArray())
    }

}


/**
 * Represents a matrix, has methods for certain matrix arithmetic
 */
export class Matrix {
    m = []
    rows = 0;
    cols = 0;

    /**
     * Make a matrix.
     * @param {number[][]} m an array of matrix rows
     */
    constructor(m) {
        this.m = m;
        this.rows = m.length;
        let t = m[0].length;

        m.forEach((arr) => { if (arr.length != t) throw "not all rows are the same length"; });

        this.cols = t;
    }

    toArray() {
        let arr = [];
        this.m.forEach((a) => { 
            let sub = []; a.forEach((b) => { sub.push(b); });
            arr.push(sub);
        })

        return arr;
    }

    toVector() {
        if (!(cols == 1 || rows == 1)) throw "matrix has more than one dimension!";
        else if (rows == 1) {
            let arr = [];
            this.m[0].forEach((a) => { arr.push(a); });
            return new Vector(arr);
        }
        else {
            let arr = [];
            this.m.forEach((a) => { arr.push(a[0]); });
            return new Vector(...arr);
        }
    }

    /**
     * Calculate the transpositon of this
     */
    T() {
        let arr = [];
        for (let i = 0; i < this.cols; i++)
            arr.push(this.colToArray(i));
        return new Matrix(arr);
    }

    /**
     * Turn a col to a 1d array
     * @param {number} c 
     */
    colToArray(c) {
        let arr = [];
        this.m.forEach((n) => { arr.push(n[c]); });
        return arr;
    }

    /**
     * Turn a row into a 1d array
     * @param {number} c 
     */
    rowToArray(c) {
        let arr = [];
        this.m[c].forEach((n) => { arr.push(n); });
        return arr;
    }

    /**
     * multiply by a matrix or a constant
     * @param {Matrix} a
     * @returns {Matrix}
     */
    mul(a) {
        //let arr = [];
        if (!(a instanceof Matrix)) {  // grab the array if is matrix
            if (Array.isArray(a)) return this.mul(new Matrix(a));  // try implicit conversion if is array
            else if (typeof a == 'number') {  // if is number do simple multiple
                // do a scalar multiply
                let res = [];
                this.m.forEach((row) => {
                    let r = [];
                    row.forEach((n) => { r.push(n*a); });
                    res.push(r);
                });
                return res;
            }
            else throw 'a must be a matrix, a number or an implicitly convertable array';
        }

        // validate if the matrices are of the correct size
        if (this.cols != a.rows) throw 'matrices must be of compatible size'

        // computer the result
        let res = [];
        for (let i = 0; i < this.rows; i++) {
            let r = [];
            for (let j = 0; j < a.cols; j++) {
                let sum = 0;
                for (let k = 0; k < this.cols; k++)
                    sum += (this.m[i][k] * a.m[k][j]);
                r.push(sum);
            }
            res.push(r);
        }

        return new Matrix(res);
    }

    /**
     * Add a matrix
     * @param {Matrix} a
     * @returns {Matrix}
     */
    add(a) {
        if (Array.isArray(a)) return this.add(new Matrix(a));
        else if (!(a instanceof Matrix)) throw 'a must be a matrix or implicitly convertable array';

        if (!(this.rows == a.rows && this.cols == a.cols)) throw 'matrices must be of the same size';

        let res = [];
        for (let i = 0; i < a.rows; i++) {
            let r = [];
            for (let j = 0; j < a.cols; j++) r.push(this.m[i][j] + a.m[i][j]);
            res.push(r);
        }

        return new Matrix(res);
    }

    /**
     * Subtract a matrix
     * @param {Matrix} a 
     */
    sub(a) {
        if (Array.isArray(a)) return this.add(new Matrix(a));
        else if (!(a instanceof Matrix)) throw 'a must be a matrix or implicitly convertable array';

        if (!(this.rows == a.rows && this.cols == a.cols)) throw 'matrices must be of the same size';

        let res = [];
        for (let i = 0; i < a.rows; i++) {
            let r = [];
            for (let j = 0; j < a.cols; j++) r.push(this.m[i][j] - a.m[i][j]);
            res.push(r);
        }

        return new Matrix(res);
    }

    /**
     * Compute the QR decomposition of this matrix.
     * For explanation see: https://en.wikipedia.org/wiki/QR_decomposition#Using_the_Gram%E2%80%93Schmidt_process
     * @returns {Matrix[]} [Q, R];
     */
    gramSchmidtQR() {

        /**
         * a useful projection
         * @param u {Vector}
         * @param a {Vector}
         * @returns {Vector}
         */
        const _proj = (u, a) => {
            //console.log('u', u, 'a', a);
            //let ua = ;
            //console.log('ua', ua);
            //let uu = ;
            //console.log('uu', uu);
            return u.mul(u.scalarProd(a)/u.scalarProd(u)); 
        };

        /** @type {Vector[]} array of colvects*/
        let A = new Array(this.cols);
        for (let i = 0; i < this.cols; i++)
            A[i] = new Vector(this.colToArray(i));
        
        //console.log('A', A);

        /** @type {Vector[]} */
        let U = new Array(A.length);
        /** @type {Vector[]} */
        let E = new Array(A.length);
        for (let i = 0; i < A.length; i++) {
            U[i] = Vector.clone(A[i]);
            //console.log(U);
            for (let j = 0; j < i; j++)
                U[i] = U[i].sub( _proj(U[j], A[i]) ); 
            E[i] = U[i].div(U[i].abs());
        }

        //console.log('U', U);
        //console.log('E', E);

        // make a matrix from the col vects in E
        let Q = [];
        for (let i = 0; i < E[0].order; i++) {
            let row = [];
            for (let j = 0; j < E.length; j++) 
                row.push(E[j].v[i]);
            Q.push(row);
        }

        //console.log('Q', Q);

        // make R
        let matQ = new Matrix(Q);
        let matR = matQ.T().mul(this)
        
        return [matQ, matR];
    }

    /**
     * Compute the determinant of this matrix
     * @returns {number}
     */
    det() {
        if (!this.isSquare()) throw 'matrices must be square'


        // QR decompose see: https://en.wikipedia.org/wiki/QR_decomposition#Using_the_Gram%E2%80%93Schmidt_process
        let R = this.gramSchmidtQR()[1];
        //let Q = QR[0];
        //let R = QR[1];

        /*
        det(A) = det(Q) * det(R)
        Q is unitrary -> det(Q) = 1
        det(A) = det(R)
        R is diagonal -> det(R) = prod(R[i][i]) // product along diagonal
        */
        let res = 1;
        for (let i = 0; i < R.cols; i++)
            res *= R.m[i][i];
        
        return res;
    }

    /**
     * Carry out gaussian elimination upon the 
     * augmented matrix [this|a].
     * @param {Matrix} a 
     * @returns {Matrix[]} the augmented matrix that is the result
     * @todo perhaps add a return that relates the det of returned triangular matrix to the det of this
     */
    gaussianElim(a) {
        // see: https://en.wikipedia.org/wiki/Gaussian_elimination#Applications

        //type checking
        if (Array.isArray(a)) return this.gaussianElim(new Matrix(a));
        else if (a instanceof Matrix == false) throw 'a must be a matrix or implicitly convertable array';
          
        // the left an right sides of the augmented matrix
        let left = Matrix.clone(this);
        let right = Matrix.clone(a);

        if (left.rows != right.rows) throw 'matrices must have the same number of rows';

        // // row ops:
        /**
         * 1. Swap the positions of two rows
         * @param {number} ra the index of first row
         * @param {number} rb the index of the second
         */
        const swap = (ra, rb) => {
            let tmp = left.m[ra];
            left.m[ra] = left.m[rb];
            left.m[rb] = tmp;
            tmp = right.m[ra];
            right.m[ra] = right.m[rb];
            right.m[rb] = tmp;
        };
        // 2. multiply a row by a non-zero scalar
        // 3. Add to one row a scalar multiple of another

        let pr = 0, pc = 0;  // pivot row, and col
        while (pr < left.rows && pc < left.cols) {

            let cnt = 1;
            while (left.m[pr][pc] == 0 && cnt < left.rows) {
                swap(pr, (pr+cnt)%left.rows);
                cnt ++;
            }
    
            for (let i = pr+1; i < left.rows; i++) {
                let f = left.m[i][pc] / left.m[pr][pc]; // work out the factor required for elimination of all x in subsequent rows
                
                left.m[i][pc] = 0;  // may improve floating point errors
                for(let j = pc+1; j < left.cols; j++)
                    left.m[i][j] -= left.m[pr][j]*f;
                for(let j = 0; j < right.cols; j++)
                    right.m[i][j] -= right.m[pr][j]*f;
            }
            pr ++;
            pc ++; 
            
        }

        return [left, right];

    }

    /**
     * Carry out Gauss-Jordan elimination upon the 
     * augmented matrix [this|a].
     * @param {Matrix} a 
     * @returns {Matrix[]} the augmented matrix that is the result
     */
    gaussJordanElim(a) {
       
        //type checking
        if (Array.isArray(a)) return this.gaussianElim(new Matrix(a));
        else if (a instanceof Matrix == false) throw 'a must be a matrix or implicitly convertable array';
          
        // the left an right sides of the augmented matrix
        let left = Matrix.clone(this);
        let right = Matrix.clone(a);

        if (left.rows != right.rows) throw 'matrices must have the same number of rows';

        /**
         * 1. Swap the positions of two rows
         * @param {number} ra the index of first row
         * @param {number} rb the index of the second
         */
        const swap = (ra, rb) => {
            let tmp = left.m[ra];
            left.m[ra] = left.m[rb];
            left.m[rb] = tmp;
            tmp = right.m[ra];
            right.m[ra] = right.m[rb];
            right.m[rb] = tmp;
        };

        let pr = 0, pc = 0;  // pivot row, and col
        while (pr < left.rows && pc < left.cols) {

            let cnt = 1;
            while (left.m[pr][pc] == 0 && cnt < left.rows) {
                swap(pr, (pr+cnt)%left.rows);
                cnt ++;
            }
    
            for (let i = 0; i < left.rows; i++) {
                if (i != pr) {
                    let f = left.m[i][pc] / left.m[pr][pc]; // work out the factor required for elimination of all x in subsequent rows
                    
                    left.m[i][pc] = 0;  // may improve floating point errors
                    for(let j = pc+1; j < left.cols; j++)
                        left.m[i][j] -= left.m[pr][j]*f;
                    for(let j = 0; j < right.cols; j++)
                        right.m[i][j] -= right.m[pr][j]*f;
                }
            }
            pr ++;
            pc ++; 
            
        }

        // now divide each row through to make left an identity matrix
        for (let i = 0; i < Math.min(left.rows, left.cols); i++) {
            for(let j = 0; j < right.cols; j++)
                right.m[i][j] /= left.m[i][i];
            left.m[i][i] = 1;
        }

        return [left, right];
    }

    /**
     * Compute the multiplicative iverse of the this matrix
     * @returns {Matrix}
     */
    inv() {
        // see: https://en.wikipedia.org/wiki/Gaussian_elimination#Applications
        if (!this.isSquare()) throw "matrix must be square to have inverse";

        let res = this.gaussJordanElim(Matrix.getI(this.rows));
        if (res[0].isI()) return res[1];
        else throw "maxtrix was uninvertable";
    }

    /**
     * Queries wheather or not the matrix is square
     * @returns {boolean}
     */
    isSquare() {
        return (this.rows == this.cols);
    }

    /**
     * Queries wheather this matrix is an identity matrix
     * @returns {boolean}
     */
    isI() {
        if (this.rows != this.cols) return false;
        
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.cols; j++) {
                //let val = (i == j)? 1 : 0;
                if (this.m[i][j] != ((i == j)? 1 : 0)) return false;
            }
        }

        return true;
    }

    /**
     * Get the identity matrix for a given size
     * @param {number} size 
     * @returns {Matrix}
     */
    static getI(size) {
        let arr = [];
        for (let i = 0; i < size; i++) {
            let row = [];
            for (let j = 0; j < size; j++) row.push((i == j)? 1 : 0);
            arr.push(row);
        }
        return new Matrix(arr);
    }
    
    /**
     * Make a clone of the passed matrix
     * @param {Matrix} a 
     */
    static clone(a) {
        let arr = [];
        a.m.forEach((r) => {
            let row = [];
            r.forEach((n) => { row.push(n); });
            arr.push(row);
        });

        return new Matrix(arr);
    }
}


