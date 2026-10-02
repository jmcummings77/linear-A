import test from "node:test";
import assert from "node:assert/strict";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { createMatrixAPI } from "./matrix.mjs";

const options = process.env.LINEAR_A_WASM_MODULE
  ? { moduleUrl: pathToFileURL(resolve(process.env.LINEAR_A_WASM_MODULE)) } : {};
const { Matrix } = await createMatrixAPI(options);
function use(items, operation) { try { operation(...items); } finally { for (const item of items) item.dispose(); } }
function near(actual, expected, tolerance = 2e-14) {
  assert.ok(Number.isFinite(actual) && Math.abs(actual - expected) < tolerance, `${actual} differs from ${expected}`);
}
function rotationInvariants(r) {
  use([r.transpose()], t => use([t.multiply(r)], product => {
    for (let row = 0; row < r.rows; row++) for (let col = 0; col < r.cols; col++) near(product.get(row,col), row === col ? 1 : 0);
  }));
  near(r.determinant(), 1);
}

test("cross products preserve the left orientation and right-hand convention", () => {
  for (const left of [1,3]) for (const right of [1,3]) {
    use([new Matrix(left,3/left,[2,3,4]),new Matrix(right,3/right,[5,6,7])], (a,b) => {
      use([a.cross(b),b.cross(a),a.cross(a)], (result,reverse,zero) => {
        assert.deepEqual([result.rows,result.cols],[left,3/left]);
        assert.deepEqual([...result.toArray()],[-3,6,-3]);
        assert.deepEqual([...reverse.toArray()],[3,-6,3]);
        assert.deepEqual([...zero.toArray()],[0,0,0]);
        result.set(0,0,99);
        assert.deepEqual([...a.toArray()],[2,3,4]);
        assert.deepEqual([...b.toArray()],[5,6,7]);
      });
    });
  }
  use([new Matrix(3,1,[1,0,0]),new Matrix(3,1,[0,1,0])], (x,y) => {
    use([x.cross(y)], z => assert.deepEqual([...z.toArray()],[0,0,1]));
  });
});

test("principal rotations actively rotate column vectors with right-handed radians", () => {
  const q = Math.PI/2;
  use([Matrix.rotation2D(q),new Matrix(2,1,[1,0])], (r,x) => use([r.multiply(x)], y => { near(y.get(0,0),0); near(y.get(1,0),1); }));
  for (const [factory,input,expected] of [
    [Matrix.rotationX,[0,1,0],[0,0,1]],
    [Matrix.rotationY,[1,0,0],[0,0,-1]],
    [Matrix.rotationZ,[1,0,0],[0,1,0]],
  ]) use([factory(q),new Matrix(3,1,input)], (r,a) => use([r.multiply(a)], result => {
    [...result.toArray()].forEach((value,i) => near(value,expected[i]));
  }));
  for (const angle of [0,q,-Math.PI,4*Math.PI,1e-8,.5,Number.MAX_VALUE]) {
    use([Matrix.rotation2D(angle)], rotationInvariants);
    for (const [index,factory] of [Matrix.rotationX,Matrix.rotationY,Matrix.rotationZ].entries()) {
      const input = [0,0,0]; input[index] = 1;
      use([new Matrix(3,1,input),factory(angle)], (axis,fixed) => use([Matrix.rotationAxisAngle(axis,angle)], general => {
        rotationInvariants(fixed); rotationInvariants(general);
        [...general.toArray()].forEach((value,i) => near(value,fixed.toArray()[i]));
      }));
    }
  }
});

test("axis normalization handles huge, tiny and subnormal scales without changing inputs", () => {
  use([new Matrix(1,3,[1,1,1])], normal => use([Matrix.rotationAxisAngle(normal,.5)], expected => {
    for (const scale of [1,1e-300,1e300,Number.MIN_VALUE]) {
      use([new Matrix(3,1,[scale,scale,scale])], axis => use([Matrix.rotationAxisAngle(axis,.5)], result => {
        rotationInvariants(result);
        [...result.toArray()].forEach((value,i) => near(value,expected.toArray()[i]));
        assert.deepEqual([...axis.toArray()],[scale,scale,scale]);
      }));
    }
  }));
  use([new Matrix(3,1,[1,1,0])], axis => use([Matrix.rotationAxisAngle(axis,1e-8)], r => near(r.get(0,1)/2.5e-17,1)));
});

test("rotations fix their axis, preserve cross products and invert when the axis reverses", () => {
  use([new Matrix(3,1,[1,2,3]),new Matrix(3,1,[2,3,4]),new Matrix(3,1,[5,6,7])], (axis,a,b) => {
    use([Matrix.rotationAxisAngle(axis,.5),axis.scale(-1),a.cross(b)], (r,negative,cross) => {
      use([r.multiply(axis),r.multiply(a),r.multiply(b),r.multiply(cross),r.transpose(),Matrix.rotationAxisAngle(negative,.5)], (fixed,ra,rb,rc,transpose,inverse) => {
        [...fixed.toArray()].forEach((value,i) => near(value,axis.toArray()[i]));
        [...inverse.toArray()].forEach((value,i) => near(value,transpose.toArray()[i]));
        use([ra.cross(rb)], actual => [...actual.toArray()].forEach((value,i) => near(value,rc.toArray()[i])));
      });
    });
  });
});

test("vector and rotation failures reject invalid shapes, overflow, angles and handles", async () => {
  use([new Matrix(3,1,[1,0,0]),new Matrix(3,3),new Matrix(1,2),new Matrix(3,1),new Matrix(3,1,[Number.MAX_VALUE,0,0]),new Matrix(3,1,[0,2,0])], (axis,square,short,zero,huge,right) => {
    assert.throws(() => axis.cross(square),/dimensions/);
    assert.throws(() => short.cross(axis),/dimensions/);
    assert.throws(() => huge.cross(right),/nonfinite/);
    assert.throws(() => Matrix.rotationAxisAngle(zero,.5),/argument/);
    assert.throws(() => Matrix.rotationAxisAngle(square,.5),/dimensions/);
    for (const angle of [NaN,Infinity,-Infinity,"0.5",true]) {
      for (const factory of [Matrix.rotation2D,Matrix.rotationX,Matrix.rotationY,Matrix.rotationZ]) assert.throws(() => factory(angle),/finite/);
      assert.throws(() => Matrix.rotationAxisAngle(axis,angle),/finite/);
    }
    assert.throws(() => Matrix.rotationAxisAngle([1,0,0],.5),/Matrix/);
    const dead = new Matrix(3,1,[1,0,0]); dead.dispose();
    assert.throws(() => dead.cross(axis),/disposed/);
    assert.throws(() => axis.cross(dead),/disposed/);
    assert.throws(() => Matrix.rotationAxisAngle(dead,.5),/disposed/);
  });
  const { Matrix: OtherMatrix } = await createMatrixAPI(options);
  use([new Matrix(3,1,[1,0,0]),new OtherMatrix(3,1,[0,1,0])], (a,b) => {
    assert.throws(() => a.cross(b),/same WebAssembly instance/);
    assert.throws(() => Matrix.rotationAxisAngle(b,.5),/same WebAssembly instance/);
  });
});
