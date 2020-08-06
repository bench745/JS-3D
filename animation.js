/**
 * @todo add sorting by dist from cam ensuring the forground is drawn in the front
 * @todo fill in the Particle class
 * @todo add some cool demos
 * @move this as a lib to its own folder
 */

import {Matrix, Vector} from './linalg.js';


/**
 * A point in 3d space
 */
class Point {
    x = null;  // x ordinate
    y = null;  // y ordinate
    z = null;  // z ordinate

    /**
     * @param {number} x 
     * @param {number} y 
     * @param {number} z 
     */
    constructor(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    /**
     * Calculate an x, y, z vector from this to p
     * @param {Point} p 
     * @returns {Vector} the vector
     */
    delta(p) {
        return new Vector([p.x-this.x, p.y-this.y, p.z-this.z]);
    }

    /**
     * @returns {Vector}
     */
    toVect() { return new Vector([this.x, this.y, this.z]); }

    /**
     * @returns {number[]}
     */
    toArray() { return [this.x, this.y, this.z]; }

    /**
     * Translate the point.
     * @param {Vector} w the vector of the translation
     */
    translate(w) {
        this.x += w.v[0];
        this.y += w.v[1];
        this.z += w.v[2];
    }

    /**
     * Unimplemented.
     * Rotate the point by phi degrees in the x y plane,
     * and by theta in the z y plane. 
     * Rotate about centre
     * @param {Point} centre 
     * @param {number} phi 
     * @param {number} theta 
     */
    rotate(centre, phi, theta) { }

}


/**
 * A camera, used to spicify the vewipoint of an observer
 */
class Camera extends Point {
    forward = new Vector([0, 0, 1]);
    up = new Vector([0, 1, 0]);
    right = new Vector([1, 0, 0]);

    fovVert = 0;  // an angle in radians, describes the angel between vision centre and verticle extremes
    fovHoriz = 0;  // an angle in radians, // and horizontal extremes
}


/**
 * A set of attributes that anything must have to be considered Drawable
 */
class Drawable {

    /** 
     * An array of Points, used as the corners of this primative 
     * has a max length of three.
     * @type {Point[]}
     */
    vertexs = [];

    /**
     * A set of [x, y, z] positions this object should be drawn at x, y and 
     * z is used for depth
     * @type {number[][]}
     */
    projected = [[],[],[]];

    /**
     * A flags to show if vertexs have changed since last projection
     * @type {boolean[]}
     */
    changed = [true, true, true];

    /**
     * Used by the Scene class to work out if this drawable needs be rendered.
     * The user needn't worry
     * @type {boolean[]}
     */
    _visible = [true, true, true];
    
    /**
     * Draw this drawable
     * @param {CanvasRenderingContext2D} ctx canvas context
     * @!param {number} width canvas width
     * @!param {number} height canvas heigth
     */
    draw(ctx) {}
}


/**
 * used to draw polygons upto the size of triangle
 */
class Primative extends Drawable {
    /**
     * The colour of the primative. Any color accepted by the 
     * canvas object will do.
     * @type {string}
     */
    colour = 'rgb(0, 0, 0)';

    /**
     * A flag used to say wheather or not the shape should be filled
     * @type {boolean}
     */
    fill = true;

    /**
     * @param {CanvasRenderingContext2D} ctx 
     */
    draw(ctx) {
        if (this.fill) ctx.fillStyle = this.colour;
        else this.strokeStyle = this.colour;

        let infront = false;
        for (let i = 0; i < this.vertexs.length; i++) infront |= (this.projected[i][2] > 0);

        if (this._visible.includes(true) || (this.fill && infront)) {
            ctx.beginPath();
            ctx.moveTo(this.projected[0][0], this.projected[0][1]);
            this.projected.forEach((coord) => {
                ctx.lineTo(coord[0], coord[1]);
            });
            ctx.lineTo(this.projected[0][0], this.projected[0][1]);
        }

        if (this.fill) ctx.fill();
        else ctx.stroke();
    }

    /**
     * Translate the shape
     * @param {Vector} v a three vector [x, y, z]
     */
    translate(v) {
        this.vertexs.forEach((vert) => {
            vert.translate(v);
        });
    }

    /**
     * UNIMPLEMENTED.
     * Rotate the shape by phi degrees in the x y plane,
     * and by theta in the z y plane. 
     * Rotate about centre
     * @param {Point} centre 
     * @param {number} phi 
     * @param {number} theta 
     */
    rotate(centre, phi, theta) { }
}


/**
 * UNIMPLIMENTED
 * Can be used to display textures, that always face the camera
 */
class Particle extends Drawable {

}




class Scene {
    cam = new Camera(0, 0.5, -3);  // a camera at (x=0, y=1, z=-1) giving a nive view of the origin

    /** used to flag which objects have moved since last frame @type {boolean[]}*/
    //moved = [];

    /** @type {Drawable[]} the set of objects to display*/
    objs = [];

    /** @type {HTMLCanvasElement} */
    cvs = null; // the canvas to draw on
    /** @type {CanvasRenderingContext2D} */
    ctx = null; // the canvas context

    constructor(canvas, fovHoriz, fovVert) {
        this.cvs = canvas;
        this.ctx = canvas.getContext('2d');

        this.cam.fovVert = fovVert;
        this.cam.fovHoriz = fovHoriz;
    }

    /**
     * Add an object to the set of visible objects
     * @param  {...Drawable} objects a collection of renderable objects
     */
    addObj(... objects) {
        this.objs.push(...objects);
        //this.moved.push(true);
    }

    /**
     * map a path from the camera to the passed 
     * @param {Point} vertex
     * @returns {number[]} units in the directions of the cameras, forward, up and right vectors
     */
    locate(vertex) {
        let v = this.cam.delta(vertex);
        let vm = v.toMatrix();

        let f = this.cam.forward.toArray();
        let u = this.cam.up.toArray();
        let r = this.cam.right.toArray();
        let cm = new Matrix([
            [f[0], u[0], r[0]],
            [f[1], u[1], r[1]],
            [f[2], u[2], r[2]],
        ]);

        let ansm = cm.gaussJordanElim(vm)[1].m;

        // could use this as return values?
        //let ansf = ansm[0][0] * this.cam.forward.abs();  // how many units in the dir of forward vect
        //let ansu = ansm[1][0] * this.cam.up.abs();  // how many units in the dir of up vects 
        //let ansr = ansm[2][0] * this.cam.right.abs();  // how many units in the dir of right vects 
        // --> to get to the this.objs[ind]
        
        return [
            ansm[0][0] * this.cam.forward.abs(), 
            ansm[1][0] * this.cam.up.abs(), 
            ansm[2][0] * this.cam.right.abs()
        ];
    }

    /**
     * project the objects vertexes onto a 2d coordinate axis
     * with -1 < x < 1 and -1 < y < 1
     */
    project() {
        for (let j = 0; j < this.objs.length; j++) {
            let item = this.objs[j];
            for (let i = 0; i < item.vertexs.length && i < 3; i++) {
                if (item.changed[i]) {  // if the drawable has moved position
                    let fur = this.locate(item.vertexs[i]);
                    let htheta = Math.atan(fur[2]/fur[0]); // verticle theta
                    let vtheta = Math.atan(fur[1]/fur[0]); // horizontal theta

                    //console.log('fur', fur, 'horiz theta', htheta, 'verticle theta', vtheta);

                    if (Math.abs(htheta) < this.cam.fovHoriz && 
                        Math.abs(vtheta) < this.cam.fovVert) item._visible[i] = true;
                    else item._visible[i] = false;

                    // this could be changed to be xy scaled by distance?
                    let w = this.cvs.width/2;
                    let h = this.cvs.height/2;
                    item.projected[i] = [
                        w + (htheta/this.cam.fovHoriz * w), 
                        h - (vtheta/this.cam.fovVert * h), 
                        fur[2]
                    ];
                    
                    //console.log('vertex visible', item._visible[i], item.projected[i]);
                }
            }
            
        }
    }

    /**
     * Render the scene 
     * @todo perhaps move the actual drawing inside the drawable, allowing for greater ability to customize
     * @todo reorder this.objs such that things in the back go in the back
     */
    render() {
        // work out projections
        this.project();

        // clear screen 
        this.ctx.clearRect(0, 0, this.cvs.width, this.cvs.height); // clear the screen

        // display projections
        this.objs.forEach( /** @param {Drawable} item */ (item) => {
            //console.log('drawable', item);
            item.draw(this.ctx);
        });

        // request next render
        //window.requestAnimationFrame(this.render);
    }


}



// the canvas and related objects
const canvas = document.getElementById('scene');  // the html obj
//const ctx = canvas.getContext('2d');  // the object used for drawing

// canvas dimensions
let width = canvas.offsetWidth;
let height = canvas.offsetHeight;

//let pxperunit = 10;
//let projmidx = width/2;
//let projmidy = height/2;

/**
 * Resize the canvas element based on the css calculated values
 */
function resize() {
    
    // redefine width and height
    width = window.innerWidth; //canvas.offsetWidth;
    height = window.innerHeight; //canvas.offsetHeight;
    
    if (height > width * (3/4)) height = width * (3/4);  // scale to meet the ratio of the fov
    //projmidx = width/2;
    //projmidy = height/2;

    // reset the canvas element's size
    canvas.width = width;
    canvas.height = height;
    //*/

}

// add the window event listener
window.addEventListener('resize', resize);
resize();  // resize the can vas to match the window


var view = new Scene(canvas, 2*Math.PI/5, 1.5*Math.PI/5);

let tri = new Primative();
tri.fill = false;
tri.vertexs = [
    new Point(-2, 0, 0),
    new Point(0, 2, 0),
    new Point(2, 0, 0)
];
view.addObj(tri);


/**
 * rotate the tri obj through its central point
 * (by moving the two extreme vertexs)
 * @param {number} theta
 */

function rotatetri(theta) {
    let centre = new Point(0, 0, 1);  // the point at the centre of the lower line
    let pta = tri.vertexs[0];
    let ptb = tri.vertexs[2];

    let cos = Math.cos(theta);
    let sin = Math.sin(theta);

    let npta = [
        cos * pta.x + sin * pta.z, 
        0, //pta.y,
        -sin * pta.x + cos * pta.z
    ];
    let nptb = [
        cos * ptb.x + sin * ptb.z, 
        0, //ptb.y,
        -sin * ptb.x + cos * ptb.z
    ];

    tri.vertexs[0] = new Point(...npta);
    tri.vertexs[2] = new Point(...nptb);
}


let hndl = setInterval(() => {
    //tri.translate(new Vector([0, 0, 0]));
    rotatetri(0.01);
    view.render();

    console.log(...tri.vertexs);
}, 10);

