
camera {
	location <-770,150,524>
	look_at <0,0,0>
	angle 20
	right x*image_width/image_height
	 }

light_source {
	<0,0,1000>
	rgb <1,1,1>
	shadowless
	}
light_source {
	<0,-1000,0>
	rgb <1,1,1>
	shadowless
	}
light_source {
	<-1000,0,0>
	rgb <1,1,1>
	shadowless
	}


background {
	color rgb <1,	1,	1>	}

difference {
	intersection {
	intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <1,1,0,0.8> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	sphere {
	<24.0276,0,13.7827>,
	36
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	box {
	<-29.2,-29.2,-29.2>,<29.2,29.2,13.0364>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-36,-36,-36>,<36,36,20.8086>
	 rotate <0,90.0002,0>
	 rotate <0,0,150>
	 translate <24.0276,0,13.7827>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-36,-36,-36>,<36,36,20.8086>
	 rotate <0,90.0002,0>
	 rotate <0,0,-150>
	 translate <24.0276,0,13.7827>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-36,-36,-36>,<36,36,19.5868>
	 rotate <0,164.839,0>
	 rotate <0,0,180>
	 translate <24.0276,0,13.7827>
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,103.186,0>
	 rotate <0,0,-1.34774>
	 translate <25.5022,-0.599988,-5.97647>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,127.862,0>
	 rotate <0,0,-102.89>
	 translate <-4.61425,-20.1633,-16.0807>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,75.8157,0>
	 rotate <0,0,-61.6765>
	 translate <12.0516,-22.3603,6.4201>
	pigment { color rgbt <1,1,0,0.8> }
	}

	cylinder {
	<2.06411,3.37825,-0.571722>,<15.068,24.6612,-4.17356>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	cylinder {
	<1.21701,-3.60138,-1.24458>,<8.88413,-26.29,-9.08544>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	cylinder {
	<3.43387,1.82299,-0.940864>,<25.0673,13.3078,-6.86832>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

difference {
	intersection {
	intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <1,1,0,0.8> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	sphere {
	<-12.0138,20.8085,13.7827>,
	36
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	box {
	<-29.2,-29.2,-29.2>,<29.2,29.2,13.0364>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-36,-36,-36>,<36,36,20.8086>
	 rotate <0,90.0002,0>
	 rotate <0,0,-30>
	 translate <-12.0138,20.8085,13.7827>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-36,-36,-36>,<36,36,20.8086>
	 rotate <0,90.0002,0>
	 rotate <0,0,-90.0002>
	 translate <-12.0138,20.8085,13.7827>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-36,-36,-36>,<36,36,25.14>
	 rotate <0,138.764,0>
	 rotate <0,0,-38.8911>
	 translate <-12.0138,20.8085,13.7827>
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,103.186,0>
	 rotate <0,0,-1.34774>
	 translate <25.5022,-0.599988,-5.97647>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,89.0348,0>
	 rotate <0,0,-177.937>
	 translate <-26.1793,-0.942829,0.441321>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,75.8157,0>
	 rotate <0,0,-61.6765>
	 translate <12.0516,-22.3603,6.4201>
	pigment { color rgbt <1,1,0,0.8> }
	}

	cylinder {
	<-3.27023,-0.142546,-2.29898>,<-23.8727,-1.04059,-16.7825>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	cylinder {
	<-2.09603,1.96014,-2.78649>,<-15.301,14.309,-20.3414>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	cylinder {
	<2.06411,3.37825,-0.571722>,<15.068,24.6612,-4.17356>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	cylinder {
	<-2.87901,1.07383,-2.56089>,<-21.0168,7.83896,-18.6945>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

difference {
	intersection {
	intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <1,1,0,0.8> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	sphere {
	<-12.0138,-20.8085,13.7827>,
	36
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	box {
	<-29.2,-29.2,-29.2>,<29.2,29.2,13.0364>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-36,-36,-36>,<36,36,20.8086>
	 rotate <0,90.0002,0>
	 rotate <0,0,30>
	 translate <-12.0138,-20.8085,13.7827>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-36,-36,-36>,<36,36,20.8086>
	 rotate <0,90.0002,0>
	 rotate <0,0,90.0002>
	 translate <-12.0138,-20.8085,13.7827>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-36,-36,-36>,<36,36,25.14>
	 rotate <0,138.764,0>
	 rotate <0,0,38.8911>
	 translate <-12.0138,-20.8085,13.7827>
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,103.186,0>
	 rotate <0,0,-1.34774>
	 translate <25.5022,-0.599988,-5.97647>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,127.862,0>
	 rotate <0,0,-102.89>
	 translate <-4.61425,-20.1633,-16.0807>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,89.0348,0>
	 rotate <0,0,-177.937>
	 translate <-26.1793,-0.942829,0.441321>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,75.8157,0>
	 rotate <0,0,-61.6765>
	 translate <12.0516,-22.3603,6.4201>
	pigment { color rgbt <1,1,0,0.8> }
	}

	cylinder {
	<-3.27023,-0.142546,-2.29898>,<-23.8727,-1.04059,-16.7825>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	cylinder {
	<1.21701,-3.60138,-1.24458>,<8.88413,-26.29,-9.08544>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	cylinder {
	<-2.09968,-2.27186,-2.53574>,<-15.3277,-16.5846,-18.5109>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

difference {
	intersection {
	intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <1,1,0,0.8> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	sphere {
	<13.7827,0,-24.0276>,
	36
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	box {
	<-29.2,-29.2,-29.2>,<29.2,29.2,13.0364>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-36,-36,-36>,<36,36,19.5868>
	 rotate <0,15.1606,0>
	 rotate <0,0,0>
	 translate <13.7827,0,-24.0276>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-36,-36,-36>,<36,36,25.14>
	 rotate <0,41.2365,0>
	 rotate <0,0,141.109>
	 translate <13.7827,0,-24.0276>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-36,-36,-36>,<36,36,25.14>
	 rotate <0,41.2365,0>
	 rotate <0,0,-141.109>
	 translate <13.7827,0,-24.0276>
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,103.186,0>
	 rotate <0,0,-1.34774>
	 translate <25.5022,-0.599988,-5.97647>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,127.862,0>
	 rotate <0,0,-102.89>
	 translate <-4.61425,-20.1633,-16.0807>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,89.0348,0>
	 rotate <0,0,-177.937>
	 translate <-26.1793,-0.942829,0.441321>
	pigment { color rgbt <1,1,0,0.8> }
	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,75.8157,0>
	 rotate <0,0,-61.6765>
	 translate <12.0516,-22.3603,6.4201>
	pigment { color rgbt <1,1,0,0.8> }
	}

	cylinder {
	<-3.27023,-0.142546,-2.29898>,<-23.8727,-1.04059,-16.7825>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	cylinder {
	<-2.09603,1.96014,-2.78649>,<-15.301,14.309,-20.3414>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	cylinder {
	<2.06411,3.37825,-0.571722>,<15.068,24.6612,-4.17356>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	cylinder {
	<1.21701,-3.60138,-1.24458>,<8.88413,-26.29,-9.08544>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	cylinder {
	<3.43387,1.82299,-0.940864>,<25.0673,13.3078,-6.86832>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	cylinder {
	<-2.87901,1.07383,-2.56089>,<-21.0168,7.83896,-18.6945>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	cylinder {
	<-2.09968,-2.27186,-2.53574>,<-15.3277,-16.5846,-18.5109>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

intersection {
	intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <0,1,1,0.8> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <0,1,1,0.8> }
	}

	}

	box {
	<-29.2,-29.2,-29.2>,<29.2,29.2,13.0364>
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,103.186,0>
	 rotate <0,0,-1.34774>
	 translate <25.5022,-0.599988,-5.97647>
	pigment { color rgbt <0,1,1,0.8> }
	}

	}

intersection {
	intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <0,1,1,0.8> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <0,1,1,0.8> }
	}

	}

	box {
	<-29.2,-29.2,-29.2>,<29.2,29.2,13.0364>
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,127.862,0>
	 rotate <0,0,-102.89>
	 translate <-4.61425,-20.1633,-16.0807>
	pigment { color rgbt <0,1,1,0.8> }
	}

	}

intersection {
	intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <0,1,1,0.8> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <0,1,1,0.8> }
	}

	}

	box {
	<-29.2,-29.2,-29.2>,<29.2,29.2,13.0364>
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,89.0348,0>
	 rotate <0,0,-177.937>
	 translate <-26.1793,-0.942829,0.441321>
	pigment { color rgbt <0,1,1,0.8> }
	}

	}

intersection {
	intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <0,1,1,0.8> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <0,1,1,0.8> }
	}

	}

	box {
	<-29.2,-29.2,-29.2>,<29.2,29.2,13.0364>
	pigment { color rgbt <1,1,0,0.8> }
	}

	}

	box {
	<-4.06,-4.06,-4.06>,<4.06,4.06,4.06>
	 rotate <0,75.8157,0>
	 rotate <0,0,-61.6765>
	 translate <12.0516,-22.3603,6.4201>
	pigment { color rgbt <0,1,1,0.8> }
	}

	}

intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <0,0,0,0.5> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

	cylinder {
	<-20.6024,-0.898041,-14.4836>,<-23.8727,-1.04059,-16.7825>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <0,0,0,0.5> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

	cylinder {
	<-13.205,12.3489,-17.5549>,<-15.301,14.309,-20.3414>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <0,0,0,0.5> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

	cylinder {
	<13.0039,21.2829,-3.60184>,<15.068,24.6612,-4.17356>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <0,0,0,0.5> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

	cylinder {
	<7.66713,-22.6887,-7.84085>,<8.88413,-26.29,-9.08544>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <0,0,0,0.5> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

	cylinder {
	<21.6334,11.4849,-5.92745>,<25.0673,13.3078,-6.86832>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <0,0,0,0.5> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

	cylinder {
	<-18.1377,6.76513,-16.1336>,<-21.0168,7.83896,-18.6945>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

intersection {
	difference {
	sphere {
	<0,0,0>,
	29.2
	pigment { color rgbt <0,0,0,0.5> }
	}

	sphere {
	<0,0,0>,
	26.2
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

	cylinder {
	<-13.228,-14.3127,-15.9752>,<-15.3277,-16.5846,-18.5109>,2.99
	pigment { color rgbt <0,0,0,0.5> }
	}

	}

  background { color rgb <0.9, 0.9, 0.9> }
sphere {
	<25.1344,-6.08962,4.19749>,
	1.25
	pigment { color rgbt <1,0,0,0.651441> }
	}

