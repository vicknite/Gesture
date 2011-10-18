function changeX(pt){
    pt.x = 10;
}

function changeX2(pt){

    pt.setX(89);
}

function write(obj){
document.write(obj +  "<br />");

}

function test(){

    var pt = new Point2D(7,5);

    write(pt.x);
    changeX(pt);
    write(pt.x);
    changeX2(pt);
    write(pt.x);
}
