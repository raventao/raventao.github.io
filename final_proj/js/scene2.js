rooms.scene2 = function() {

lib3D2();

description = `
               
        <h1>Final project</h1>
        <h3>Cub bears dancing around a Christmas tree</h3>
        
        <br>
        <br>
        Ruixiang Han
        <br>
        rh2981@nyu.edu                
        `;

code = {
'init':`


   // DEFINE MATERIALS TO BE RENDERED VIA PHONG REFLECTANCE MODEL
S.black = [0,0,0,0,      0,0,0,0,    2,2,2,20,0,0,0,0]; // BLACK 
S.canopyColor = [.0,.06,.045,0, .03,.074,.04,.1, .4,.3,.1,13, 0,0,0,0];// canopy color
S.copper = [.15,.05,.025,0, .3,.1,.05,0, .6,.2,.1,15, 0,0,0,0];// COPPER
S.bearColor = [.18,.08,.06,0, .12,.08,.05,.1, .2,.3,.1,20, 0,0,0,0];// bear color
S.trunkColor = [.10,.03,.05,0, .12,.08,.05,.1, .2,.3,.1,20, 0,0,0,0];// bear color
S.gold =[.25,.15,.025,0, .5,.3,.05,0, 1,.6,.1,6,  0,0,0,0]; // GOLD
S.red_plastic = [.25,0,0,0,      .5,0,0,0,    2,2,2,20,0,0,0,0]; // PLASTIC
S.white_plastic = [1,1,1,0,      .25,0,0,0,    2,2,2,20,0,0,0,0]; // WHITE PLASTIC
S.black_plastic = [0,0,0,0,      0,0,0.1,0,    2,2,2,20,0,0,0,0]; // BLACK PLASTIC
S.lead = [.05,.05,.05,0,  .1,.1,.1,0,  1,1,1,5,    0,0,0,0]; // LEAD
S.silver = [.1,.1,.1,0,     .1,.1,.1,0,  1,1,1,5,    0,0,0,0]// SILVER;

   // A SQUARE IS A TRIANGLE MESH WITH JUST TWO TRIANGLES

   S.squareMesh = [ -1, 1, 0,  0,0,1,  0,1,
                     1, 1, 0,  0,0,1,  1,1,
                -1,-1, 0,  0,0,1,  0,0,
                 1,-1, 0,  0,0,1,  1,0 ];

   // OCTAHEDRON MESH
   S.octahedronMesh = [
      0,0,1, 1,1,1, 0,0,
      1,0,0, 1,1,1, 0,0,
      0,1,0, 1,1,1, 0,0,

      0,0,1, -1,1,1, 0,0,
      -1,0,0, -1,1,1, 0,0,
      0,1,0, -1,1,1, 0,0,

      0,0,1, -1,-1,1, 0,0,
      -1,0,0, -1,-1,1, 0,0,
      0,-1,0, -1,-1,1, 0,0,

      0,0,1, 1,-1,1, 0,0,
      1,0,0, 1,-1,1, 0,0,
      0,-1,0, 1,-1,1, 0,0,

      0,0,-1, 1,1,-1, 0,0,
      1,0,0, 1,1,-1, 0,0,
      0,1,0, 1,1,-1, 0,0,

      0,0,-1, -1,1,-1, 0,0,
      -1,0,0, -1,1,-1, 0,0,
      0,1,0, -1,1,-1, 0,0,

      0,0,-1, -1,-1,-1, 0,0,
      -1,0,0, -1,-1,-1, 0,0,
      0,-1,0, -1,-1,-1, 0,0,

      0,0,-1, 1,-1,-1, 0,0,
      1,0,0, 1,-1,-1, 0,0,
      0,-1,0, 1,-1,-1, 0,0,
   ];

   S.octahedronMesh.isTriangles = true;


   // GLUE TOGETHER TWO MESHES TO CREATE A SINGLE MESH

   let glueMeshes = (a,b) => {
      let mesh = a.slice();
      mesh.push(a.slice(a.length - S.VERTEX_SIZE, a.length));
      mesh.push(b.slice(0, S.VERTEX_SIZE));
      mesh.push(b);
      return mesh.flat();
   }

   // GIVEN A FUNCTION THAT MAPS (u,v) TO point AND normal,
   // AND GIVEN A MESH RESOLUTION, CREATE A PARAMETRIC MESH

   let uvMesh = (f,nu,nv) => {
      let mesh = [];
      for (let iv = 0 ; iv < nv ; iv++) {
         let v = iv / nv;
  let strip = [];
         for (let iu = 0 ; iu <= nu ; iu++) {
     let u = iu / nu;
     strip = strip.concat(f(u,v));
     strip = strip.concat(f(u,v+1/nv));
  }
  mesh = glueMeshes(mesh, strip);
      }
      return mesh;
   }

   // CREATE A UNIT SPHERE PARAMETRIC MESH

   S.sphereMesh = uvMesh((u,v) => {
      let theta = 2 * Math.PI * u;
      let phi = Math.PI * v - Math.PI/2;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      let cv = Math.cos(phi);
      let sv = Math.sin(phi);
      return [cu * cv, su * cv, sv,
              cu * cv, su * cv, sv,
       u, v];
   }, 20, 10);

   // CREATE A UNIT CONE WITH CAP PARAMETRIC MESH
   S.coneMesh = uvMesh((u,v) => {
      let theta = 2 * Math.PI * u;
      let phi = Math.PI * v - Math.PI/2;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      return [cu * v, su * v, 2 * v - 1,
              cu    , su    , -u / 2   ,
              u     , v];
   }, 20, 10);
   
   // CREATE A UNIT CAP PARAMETRIC MESH
   
   S.diskMesh = uvMesh((u, v) => {
      let theta = 2 * Math.PI * u;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      return [v * cu, v * su, 0, 
            0, 0, 1, 
            u, v];
   }, 20, 2);

   // CREATE A UNIT OPEN TUBE PARAMETRIC MESH

   S.tubeMesh = uvMesh((u,v) => {
      let theta = 2 * Math.PI * u;
      let phi   = 2 * Math.PI * v;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      return [cu, su, 2 * v - 1,   cu, su, 0,   u, v];
   }, 20, 2);

   // TRANSFORM A MESH BY A MATRIX ON THE CPU

   let transformMesh = (mesh, matrix) => {
      let result = [];
      let IMT = matrixTranspose(matrixInverse(matrix));
      for (let n = 0 ; n < mesh.length ; n += S.VERTEX_SIZE) {
         let V = mesh.slice(n, n + S.VERTEX_SIZE);
  let P  = V.slice(0, 3);
  let N  = V.slice(3, 6);
  let UV = V.slice(6, 8);
  P = matrixTransform(matrix, [P[0], P[1], P[2], 1]);
  N = matrixTransform(IMT,    [N[0], N[1], N[2], 0]);
         result.push(P[0],P[1],P[2], N[0],N[1],N[2], UV);
      }
      return result.flat();
   }

   // A CYLINDER MESH IS A TUBE WITH TWO DISK END-CAPS GLUED TOGETHER

   let end0 = transformMesh(S.diskMesh, matrixTranslate([0,0,1]));
   let end1 = transformMesh(end0      , matrixRotx(Math.PI));
   S.cylinderMesh = glueMeshes(S.tubeMesh, glueMeshes(end0, end1));


   // GLUE CAP TO CONE
   S.cone = glueMeshes(S.coneMesh, transformMesh(S.diskMesh, matrixTranslate([0,0,1])));


   // A CUBE MESH IS SIX TRANSFORMED SQUARE MESHES GLUED TOGETHER

   let face0 = transformMesh(S.squareMesh, matrixTranslate([0,0,1]));
   let face1 = transformMesh(face0,        matrixRotx( Math.PI/2));
   let face2 = transformMesh(face0,        matrixRotx( Math.PI  ));
   let face3 = transformMesh(face0,        matrixRotx(-Math.PI/2));
   let face4 = transformMesh(face0,        matrixRoty(-Math.PI/2));
   let face5 = transformMesh(face0,        matrixRoty( Math.PI/2));
   S.cubeMesh = glueMeshes(face0,
                glueMeshes(face1,
                glueMeshes(face2,
                glueMeshes(face3,
                glueMeshes(face4,
                 face5)))));

   // DRAW A SINGLE MESH WITH MATERIAL PROPERTIES

   S.drawMesh = (mesh, matrix, material) => {
      let gl = S.gl;
      S.setUniform('Matrix4fv', 'uMatrix', false, matrix);
      S.setUniform('Matrix4fv', 'uInvMatrix', false, matrixInverse(matrix));
      S.setUniform('Matrix4fv', 'uMaterial', false, material);
      S.gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(mesh), gl.STATIC_DRAW);
      // allow to choose between triangle and triangel strips
      S.gl.drawArrays(mesh.isTriangles ? S.gl.TRIANGLES
                              : S.gl.TRIANGLE_STRIP, 0, mesh.length / S.VERTEX_SIZE);
   }

`,
fragment: `
S.setFragmentShader(\`
   // DECLARE CONSTANTS, UNIFORMS, VARYING VARIABLES
   const int nL = \` + S.nL + \`;
   uniform vec3 uBgColor;
   uniform mat4 uMaterial;
   uniform float uTime;
   uniform vec3 uLd[nL];
   uniform vec3 uLc[nL];

   varying vec3 vPos, vNor;

   // phong shading
   // params: P - intersection point, N - surface normal, phong - material
   vec3 shadeSurface(vec3 P, vec3 N, mat4 phong){
      // extract phong params
      vec3 ambient = phong[0].rgb;
      vec3 diffuse = phong[1].rgb;
      vec3 specular = phong[2].rgb;
      float p = phong[2].a; // shininess factor

      // compute init color and approximate vector to eye
      vec3 c = mix(ambient, uBgColor, .3);//mix in ambient light for all pixels
      vec3 E = vec3(0.,0.,1.);

      // loop through light sources
      for (int l = 0; l < nL; l ++){
         vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
         c += uLc[l] * (diffuse * max(0., dot(N, uLd[l]))
                               + specular * pow(max(0., dot(R, E)),p));
      }
      //add spotty texture to copper, plastic and lead
      //if (ambient.g <= .05) c *= 1. + .5 * noise(30.*N);
      return c;
   }


   void main() {
      //float c = .2 + .8 * max(0.,dot(vNor,vec3(.57)));
      vec3 c = shadeSurface(vPos, vNor, uMaterial);
      gl_FragColor = vec4(sqrt(c), 1.);
   }
\`);
`,
vertex: `
S.setVertexShader(\`

   attribute vec3 aPos, aNor;
   varying   vec3 vPos, vNor;
   uniform   mat4 uMatrix, uInvMatrix, uProject, uMaterial;

   void main() {
      vec4 pos = uProject * uMatrix * vec4(aPos, 1.);
      vec4 nor = vec4(aNor, 0.) * uInvMatrix;
      vPos = pos.xyz;
      vNor = normalize(nor.xyz);
      gl_Position = pos * vec4(1.,1.,-.01,1.);
   }
\`)
`,
render: `

   // USEFUL VECTOR FUNCTIONS

   let add = (a,b) => [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ];
   let dot = (a,b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
   let norm = v => Math.sqrt(dot(v,v));
   let normalize = v => { let s = norm(v); return [ v[0]/s, v[1]/s, v[2]/s ]; }
   let scale = (v,s) => [ s * v[0], s * v[1], s * v[2] ];
   let subtract = (a,b) => [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ];

   // SEND LIGHT SOURCE DATA TO GPU

   let ldData = [ normalize([1,1,1]),
                  normalize([-1,-1,-1]) ];
   S.setUniform('3fv', 'uLd', ldData.flat());
   S.setUniform('3fv', 'uLc', [ 1,1,1, 1,1,1 ]);

   // DEFINE NUMBER OF LIGHTS FOR GPU

   S.nL = ldData.length;



   // SET THE PROJECTION MATRIX BASED ON CAMERA FOCAL LENGTH

   let fl = 5.0;
   S.setUniform('Matrix4fv', 'uProject', false,
      [1,0,0,0, 0,1,0,0, 0,0,1,-1/fl, 0,0,0,1]);

   let m = new Matrix();
   // tree
   m.roty(Math.PI/4*3);
   m.roty(2*Math.cos(1/4*time));
   m.save();
      m.translate([-.3,-.35,0]);

      //trunk
      m.save();
         m.scale([.07,.35,.07]);
         m.rotx(Math.PI/2);
         S.drawMesh(S.cylinderMesh, m.get(), S.trunkColor);
      m.restore();

      //canopy
      m.save();
         for (let i = 0; i<3; i++){
            m.save();
               m.translate([0,.35+i*.14,0]);
               m.scale([.3-i*.05,.3-i*.05,.3-i*.05]);
               m.rotx(Math.PI/2);
               S.drawMesh(S.coneMesh, m.get(), S.canopyColor);
            m.restore();
         }
         // dec ball
         m.scale([.05,.05,.05]);
         //1
         m.save();
            m.translate([-4.5,3.5,-3]);
            S.drawMesh(S.sphereMesh, m.get(), S.silver);
         m.restore();
         //2
         m.save();
            m.translate([4,6,.8]);
            S.drawMesh(S.sphereMesh, m.get(), S.gold);
         m.restore();
         // 3
         m.save();
            m.translate([3,9,-2.5]);
            S.drawMesh(S.sphereMesh, m.get(), S.red_plastic);
         m.restore();
         //4
         m.save();
            m.translate([-3,3.5,4.5]);
            S.drawMesh(S.sphereMesh, m.get(), S.red_plastic);
         m.restore();
         //5
         m.save();
            m.translate([-4,7,1]);
            S.drawMesh(S.sphereMesh, m.get(), S.gold);
         m.restore();
      m.restore();
   m.restore();

   // bear
   // m.roty(Math.PI/2);
   // m.roty(Math.sin(time));
   for(let i = 0;i<=1; i++){
      m.save();
         m.roty(i*Math.PI*3/2)
         m.translate([.5+i*.4,0,i*.3]);
         m.roty(i*Math.PI/4);
         m.save();
            m.translate([.5,0.003*Math.cos(time*10)+.3,0]);
            m.rotx(Math.cos(time)/20);
         m.restore();

         m.save();
         
            //head
            
               //m.translate([0,-.25,0]);
               m.rotz(Math.cos(time)/10);
               m.rotx(Math.sin(time)/10)
                  m.save();
                     m.scale([.2,.2,.2]);
                     S.drawMesh(S.sphereMesh, m.get(), S.bearColor);
                  m.restore();
                  //nose
                  m.save();
                     m.translate([-.18,.05,0]);
                     m.scale([.06,.06,.06]);
                     S.drawMesh(S.sphereMesh, m.get(), S.copper);
                  m.restore();
                  //small nose
                  m.save();
                     m.translate([-.22,.09,0]);
                     m.scale([.01,.01,.01]);
                     S.drawMesh(S.sphereMesh, m.get(), S.black);
                  m.restore();
            

            //eye

               for(let i = -1;i <=1; i+=2){
                  m.save();
                  m.translate([-.11,.10,i*(.11)]);
                  m.scale([.04,.045,.04]);
                  S.drawMesh(S.sphereMesh, m.get(), S.white_plastic);
                  m.restore();
               }
               

            //iris

               for(let i = -1;i <=1; i+=2){
                  m.save();
                     m.translate([-.13,0.11,i*(-0.11)]);
                     m.scale([.03,0.02,.02]);
                     S.drawMesh(S.sphereMesh, m.get(), S.black_plastic);
                  m.restore();
               }

            //ears
               m.save();
                  m.translate([.05,.25,.1]);
                  for(let i = -1;i <=1; i+=2){
                     m.save();
                        m.translate([.05,-0.125,-.11+i*0.165]);
                        m.rotz(-0.2);
                        m.scale([.05,.05,.05]);
                        S.drawMesh(S.sphereMesh, m.get(), S.bearColor);
                     m.restore();
                  }
               m.restore();

            //arms
                  m.save();
                  m.translate([0,.5,0]);
                  for(let i =0;i<=1;i++){
                        m.save();
                        m.translate([-.15,-.65,0-.4*i]);
                        i == 0? m.rotz(Math.sin(time+i)): m.rotz(Math.cos(time));
                        m.translate([-.15,-.05,.2]);
                        m.scale([.2,.065,.06]);
                        S.drawMesh(S.sphereMesh, m.get(), S.bearColor);
                        m.restore();
                     }
                  m.restore();

            // lower body
                  m.save();
                     m.translate([-.05,-.35,0])

                     //torso
                     m.save();
                        m.rotz(-.1)
                        m.scale([.2,.26,.23])
                        S.drawMesh(S.sphereMesh, m.get(), S.bearColor);
                        m.save();
                         m.translate([-.1,-.45,0])
                           m.scale([1,.2,.2])
                           S.drawMesh(S.sphereMesh, m.get(), S.bearColor);
                        m.restore();
                        //belly
                        m.save();
                         m.translate([-.1,-.3,0])
                           m.scale([1,.6,.7])
                           S.drawMesh(S.sphereMesh, m.get(), S.copper);
                        m.restore();
                     m.restore();
                     // feet
                     m.save();

                        m.rotx(Math.PI/2);
                        //m.roty(-Math.PI/2);
                        //m.roty(1/2*Math.cos(time));
                        for(let i = -1;i <=1; i+=2){
                           m.save()
                           m.translate([-.03,.18*i,.23]);
                           m.roty(-Math.PI/2);
                           m.scale([.09,.07,.05])
                           S.drawMesh(S.cylinderMesh, m.get(), S.bearColor);
                           m.restore();
                        }

                        // claw
                        for(let i = -1;i <=1; i+=2){

                           //the big padding
                           m.save()
                              m.translate([-.04,.18*i,.24]);
                              m.roty(-Math.PI/2);
                              m.scale([.04,.04,.055])
                              S.drawMesh(S.cylinderMesh, m.get(), S.copper);

                           m.restore();

                           // three small ones
                           for(let j = -1;j <=1; j+=1){
                              m.save()
                                 m.translate([-0.045,.18*i+j*0.03,.17+0.015*Math.abs(j)])
                                 m.roty(-Math.PI/2);
                                 m.scale([.015,.015,.055])
                                 S.drawMesh(S.cylinderMesh, m.get(), S.copper);
                              m.restore();
                           }
                        }
                     m.restore();

                     //tail
                     m.save()
                        m.translate([.18,-.15,0])
                        m.scale([.05,.05,.05])
                        S.drawMesh(S.sphereMesh, m.get(), S.bearColor);
                     m.restore();

                  m.restore();
            m.restore();
         m.restore;
      }

`,
events: `
   ;
`
};

}