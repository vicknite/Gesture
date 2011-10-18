

var Centroid = function(x,y){
    this.x = x;
    this.y = y;
}

//------------------Point2D--------------
var Point2D = function (x,y) {
    this.x = x;
    this.y = y;
}



Point2D.prototype = {
    setPoint:function(x,y){
		 this.x = x;
		 this.y = y;
	     },
    getX:function(){return this.x;},
    getY:function(){return this.y;},
    setX:function(x){this.x = x},
    setY:function(y){this.y = y}

}

//-----------------Template--------------
var Template = function (id, coll){
    this. id = id; //string as id
    this.collection = coll;
}
//----------------Collection---------------

var Collection =function (){
    this.collection ={};
    this.count = 0;
}

Collection.prototype = {
    add:function(key,item){
	    if(this.collection[key]!= undefined){return undefined;}
	    this.collection[key] = item;
	    return this.count++;
	},
    remove:function(key){
	       if(this.collection[key ]== undefined){return undefined;}
	       delete this.collection[key];
	       return this.count--;
	   },
    size:function(){return this.count},
    indexOf:function(item){},
    getItemAt:function(key){return this.collection[key];},
    clone:function(){
	      var coll = new Collection();
	      for(key in this.collection){

		  coll[key] = this.collection[key];
	      }

	      return coll;
	  },
    clear:function(){
	      for(key in this.collection){
		  delete this.collection[key];
	      }
	      this.count = 0;
	  },
    set:function(coll){
	    if (coll instanceof Collection){
		this.collection = coll.collection;
	    }
	},
}

//-----------------Pattern-----------------
var Pattern = function(template,segments){
    this.template = template;
    this.segments = segments;//list of list of Point2D
}
//-----------------Rect----------------
var Rect = function(x,y,width,height){
    this.x =x;
    this.y=y;
    this.height=height;
    this.width=width;
}

//----------------Incremental Result--------------
var IncrementalResult = function(pattern,  prob,  indexOfMostLikelySegment){
    this.pattern = pattern;
    this.prob = prob;
    this.indexOfMostLikelySegment = indexOfMostLikelySegment;
}


//-------------------Iterator-------------------------
var Iterator= function(collection){
    this.collection = collection.collection;
    this.count = collection.size();
    this.index = 0;
}

Iterator.prototype = {

    next:function(){
	     var element = this.collection[index];
	     index++;
	     return element;
	 },
    hasNext:function(){
		if(this.collection[index+1] == undefined){

		    return true;
		}
		return false;
	    },
    remove:function(){
	       if(this.collection[index] != undefined){

		   delete this.collection[index];
	       }
	       return undefined;
	   }
}

//------------------
var Result = function(template,prob,pts){
    this.template = template;	//class template
    this.prob = prob; //double prob
    this.pts = pts; //list of points
}

Result.prototype.compareTo = function(r){
    if (this.prob === r.prob) {
	return 0;
    }
    else if (prob < r.prob) {
	return 1;
    }
    else {
	return -1;
    }
}

//--------------Gesture---------------------
var Gesture = function(templates,samplePointDistance){

    this.DEFAULT_E_SIGMA = 200.0;
    this.DEFAULT_BETA = 400.0;
    this.DEFAULT_LAMBDA = 0.4;
    this.DEFAULT_KAPPA = 1.0;
    this.MAX_RESAMPLING_PTS = 1000;
    this.normalizedSpace = new Rect(0, 0, 1000, 1000);
    this.samplePointDistance =samplePointDistance ;
    this.patterns = new Collection();
    this.centroid = new Centroid(0,0);
}

Gesture.prototype = {

    deepCopy:function(pts){
		 arr = new Collection();

		 for(key in pts.collection){
		     var point = pts.getItemAt(key);
		     arr.add(key,new Point2D(point.x,point.y));
		 }

		 return arr;

	     },
    distance:function (pt1,pt2){

		 return this.distance(pt1.x,pt1.y,pt2.x,pt2.y);
	     },
    distance:function (x1,y1,x2,y2){

		 if ((x2 -= x1) < 0) {
		     x2 = -x2;
		 }
		 if ((y2 -= y1) < 0) {
		     y2 = -y2;
		 }
		 return (x2 + y2 - (((x2 > y2) ? y2 : x2) >> 1) );
	     },
    getBoundingBox:function(pts){

		       var minX,maxX,minY,maxY;
		       minX = 999999999; maxX = -999999999;
		       minY = 999999999; maxY = -999999999;

		       var i,length;
		       length = pts.size();

		       for(i = 0;i < size;i++ ){
			   var pt = pts.getItemAt(i);

			   if (pt.x < minX) {
			       minX = pt.x;
			   }
			   if (pt.x > maxX) {
			       maxX = pt.x;
			   }
			   if (pt.y < minY) {
			       minY = pt.y;
			   }
			   if (pt.y > maxY) {
			       maxY = pt.y;
			   }
		       }
		       return new Rect(minX, minY, (maxX - minX), (maxY - minY));
		   },
    getCentroid:function(pts){
		    var totalMass = pts.size();
		    var xIntegral = 0.0;
		    var yIntegral = 0.0;

		    var i,length;
		    length = pts.size();
		    for(i = 0;i < size;i++ ){
			var pt = pts.getItemAt(i);
			xIntegral += pt.x;
			yIntegral += pt.y;
		    }

		    return new Centroid(xIntegral / totalMass, yIntegral / totalMass);

		},
    getEuclideanDistance:function(pt1, pt2){

			    return Math.sqrt(this.getSquaredEuclidenDistance(pt1,pt2));
			},

    getEuclidianDistance:function (pts1, pts2){

			    if (pts1.size() != pts2.size()) {
				console.log("lists must be of equal lengths, cf. " + pts1.size() + " with " + pts2.size());
			    }
			    var n = pts1.size();
			    var td = 0;
			    for (var i = 0; i < n; i++) {
				td+= this.getEuclideanDistance(pts1.getItemAt(i), pts2.getItemAt(i));
			    }
			    return td / n;
			},
    getIncrementalResult:function(unkpts, pattern,  beta,  lambda, e_sigma){
			     var segments = pattern.segments;
			     var maxProb = 0.0;
			     var maxIndex = -1;
			     var size = segments.size();
			     for(var i = 0;i< size;i++){
				 var pts =  segments.getItemAt(i);
				 var samplingPtCount = pts.size();
				 var unkResampledPts = this.resample(unkPts, samplingPtCount);
				 var prob = this.getLikelihoodOfMatch(unkResampledPts, pts, e_sigma, e_sigma/beta, lambda);
				 if (prob > maxProb) {
				     maxProb = prob;
				     maxIndex = i;
				 }
			     }
			     return new IncrementalResult(pattern, maxProb, maxIndex);
			 },
    getIncrementalResults:function(input, beta, lambda, kappa, e_sigma){
			      var results = new Collection();
			      var unkpts = this.deepCopy(input);
			      this.normalize(unkpts);
			      for(key in this.patterns.collection){

				  var result = this.getIncrementalResult(unkPts, pattern, beta, lambda, e_sigma);
				  var lastSegmentPts = pattern.segments.getItemAt(pattern.segments.size()-1);
				  var completeProb = this.getLikelihoodOfMatch(this.resample(unkPts, lastSegmentPts.size()), lastSegmentPts, e_sigma, e_sigma/beta, lambda);
				  var x = 1 - completeProb;
				  result.prob *= (1 + kappa*Math.exp(-x*x));
				  results.add(key,result);
			      }
			      this.marginalizeIncrementalResults(results);
			      return results;
			  },
    getLikelihoodOfMatch:function(pts1,pts2,eSigma,aSigma,lambda){

			     try{
				 if (eSigma === 0 || eSigma < 0) {
				     throw 1;
				 }
				 if (aSigma === 0 || eSigma < 0) {
				     throw 2;
				 }
				 if (lambda < 0 || lambda > 1) {
				     throw 3;
				 }
				 var x_e = this.getEuclidianDistance(pts1, pts2);
				 var x_a = this.getTurningAngleDistance(pts1, pts2);


			     }catch(er){

				 switch(er){

				 case 1:
				     console.log("eSigma must be positive");
				     break;
				 case 2:
				     console.log("aSigma must be positive");
				     break;
				 case 3:
				     console.log("lambda must be in the range between zero and one")
				     break;
				 default:
				     break;
				 }
			     }

			     return Math.exp(- (x_e * x_e / (eSigma * eSigma) * lambda + x_a * x_a / (aSigma * aSigma) * (1 - lambda)));
			 },

    getResult:function(incrementalResults){
		  var results = new Collection();
		  for (key in incrementalResults.collection){
		      var ir = incrementalResults.getItemAt(key);
		      var r  = new Result(ir.pattern.template, ir.prob, ir.pattern.segments.get(ir.indexOfMostLikelySegment));
		      results.add(key,r);
		  }
		  return results;
	      },
    getSquaredEuclidenDistance:function (pt1, pt2){

				   return (pt1.x - pt2.x) * (pt1.x - pt2.x) + (pt1.y - pt2.y) * (pt1.y - pt2.y);
			       },
    getSegmentPoints:function (pts,n,length,buffer){
			 var i, m;
			 var x1, y1, x2, y2, ps;
			 var rest, currentLen;

			 m = n * 2;
			 rest = 0.0;
			 x1 = pts[0];
			 y1 = pts[1];
			 for (i = 2; i < m; i += 2) {
			     x2 = pts[i];
			     y2 = pts[i + 1];
			     currentLen = distance(x1, y1, x2, y2);
			     currentLen += rest;
			     rest = 0.0;
			     ps = (currentLen / length);
			     if (ps == 0) {
				 rest += currentLen;
			     }
			     else {
				 rest += currentLen - (ps * length);
			     }
			     if (i == 2 && ps == 0) {
				 ps = 1;
			     }
			     buffer[(i / 2) - 1] = ps;
			     x1 = x2;
			     y1 = y2;
			 }
			 return rest;

		     },
    getSpatialLength:function (pts){

			 var len = 0.0;

			 var it = Iterator(pts);

			 if(it.hasNext()){
			     var pt0 = it.next();
			     while(it.hasNext()){
				 var pt1 = it.next();
				 len += distance(pt0,pt1);
				 pt0 = pt1;
			     }
			 }
			 return len;
		     },
    getTurningAngleDistance:function (ptA1,ptA2,ptB1,ptB2){
				var len_a =this.getEuclideanDistance(ptA1,ptA2) ;
				var len_b =this.getEuclideanDistance(ptB1,ptB2) ;

				if(len_a === 0 || len_b === 0){return 0.0;}
				else{

				    var cos = (((ptA1.x - ptA2.x) * (ptB1.x - ptB2.x) + (ptA1.y - ptA2.y)*(ptB1.y - ptB2.y) ) / (len_a * len_b));
				    if(Math.abs(cos)>1.0){
					return 0.0;
				    }else{
					return Math.acos(cos);
				    }
				}
			    },
    getTurningAngleDistance:function(pts1, pts2){
				if(pts1.size() != pts2.size()){ console.log("must be equivalent size"); }
				var n = pts1.size();
				var td = 0;

				for (var i = 0; i < n - 1; i++) {
				    td+= Math.abs(getTurningAngleDistance(pts1.getItemAt(i), pts1.getItemAt(i + 1), pts2.getItemAt(i), pts2.getItemAt(i + 1)));
				}
				if ( isNaN(td)) {
				    return 0.0;
				}
				return td / (n - 1);

			    },
    marginalizeIncrementalResults:function(results){
				      var totalMass = 0.0;
				      for (key in results.collection) {
					  var r = getItemAt(key);
					  totalMass+= r.prob;
				      }

				      for (key in results.collection) {
					  var i = getItemAt(key);
					  i.prob/= totalMass;
				      }
				  },
    normalize: function(pts){
		   this.scaleTo(pts,this.normalizedSpace);
		   var c = this.getCentroid(pts);
		   this.translate(pts, -c.x, -c.y)
	       },
    recognize:function( input, beta, lambda, kappa, e_sigma){
		  if (input.size < 2){}
		  var incResults = this.getIncrementalResults(input, beta, lambda, kappa, e_sigma);
		  results = getResults(incResults);
		  this.sortResult(results,function(a,b){return a.prob - b.prob;})
		  return results;
	      },
    resample:function (points,numTargetPoints) {
		 var r = new Collection();
		 var inArray = toArray(points);
		 var outArray = new Array();

		 resample(inArray, outArray, points.size(), numTargetPoints);
		 for (var i = 0, n = outArray.length; i < n; i+= 2) {
		     r.add(i,new Point2D(outArray[i], outArray[i + 1]));
		 }
		 return r;
	     },
    resample:function(template,buffer,n,numTargetPoints) {

		 var segment_buf = Array(MAX_RESAMPLING_PTS);

		 var l, segmentLen, horizRest, verticRest, dx, dy;
		 var x1, y1, x2, y2;
		 var i, m, a, segmentPoints, j, maxOutputs, end;

		 m = n * 2;
		 l = getSpatialLength(template, n);
		 segmentLen = l / (numTargetPoints - 1);
		 getSegmentPoints(template, n, segmentLen, segment_buf);
		 horizRest = 0.0;
		 verticRest = 0.0;
		 x1 = template[0];
		 y1 = template[1];
		 a = 0;
		 maxOutputs = numTargetPoints * 2;
		 for (i = 2; i < m; i += 2) {
		     x2 = template[i];
		     y2 = template[i + 1];
		     segmentPoints = segment_buf[(i / 2) - 1];
		     dx = -1.0;
		     dy = -1.0;
		     if (segmentPoints - 1 <= 0) {
			 dx = 0.0;
			 dy = 0.0;
		     }
		     else {
			 dx = (x2 - x1) / segmentPoints;
			 dy = (y2 - y1) / segmentPoints;
		     }
		     if (segmentPoints > 0) {
			 for (j = 0; j < segmentPoints; j++) {
			     if (j == 0) {
				 if (a < maxOutputs) {
				     buffer[a] = x1 + horizRest;
				     buffer[a + 1] = y1 + verticRest;
				     horizRest = 0.0;
				     verticRest = 0.0;
				     a += 2;
				 }
			     }
			     else {
				 if (a < maxOutputs) {
				     buffer[a] = x1 + j * dx;
				     buffer[a + 1] = y1 + j * dy;
				     a += 2;
				 }
			     }
			 }
		     }
		     x1 = x2;
		     y1 = y2;
		 }
		 end = (numTargetPoints * 2) - 2;
		 if (a < end) {
		     for (i = a; i < end; i += 2) {
			 buffer[i] = (buffer[i - 2] + template[m - 2]) / 2;
			 buffer[i + 1] = (buffer[i - 1] + template[m - 1]) / 2;
		     }
		 }
		 buffer[maxOutputs - 2] = template[m - 2];
		 buffer[maxOutputs - 1] = template[m - 1];
	     },
    scaleTo:function(pts,targetBounds){
		var bounds = this.getBoundingBox(pts);
		var a1 = targetBounds.width;
		var a2 = targetBounds.height
		var b1 = bounds.width;
		var b2 = bounds.height;
		var scale = Math.sqrt(a1 * a1 + a2 * a2) / Math.sqrt(b1 * b1 + b2 * b2);
		this.scale(pts, scale, scale, bounds.x, bounds.y);
	    },
    scale:function(pts,sx,sy,originX,originY){
	      this.translate(pts, -originX, -originY);
	      this.scale(pts, sx, sy);
	      this.translate(pts, originX, originY);
	    },
    scale:function(pts,sx,sy){
	      var i,length;
	      length = pts.size();
	      for(i = 0;i < size;i++ ){
		  var pt = pts.getItemAt(i);
		  pt.x *= sx;
		  pt.y *= sy;
	      }
	    },

    sortResult:function(results,sortFunction){
		   results.collection.sort(sortFunction);
	       },
    translate:function(pts,dx,dy){

		  var i,length;
		  length = pts.size();
		  for(i = 0;i < size;i++ ){
		      var pt = pts.getItemAt(i);
		      pt.x += Math.floor(dx);
		      pt.y += Math.floor(dy);
		  }

	      },
    toArray:function (points){
		var out = new Array();
		for(key in points.collection){

		    out.push(points.collection[key].x);
		    out.push(points.collection[key].y);
		}
		return out;
	    }
}


