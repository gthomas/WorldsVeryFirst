parcelRequire=function(e,r,t,n){var i,o="function"==typeof parcelRequire&&parcelRequire,u="function"==typeof require&&require;function f(t,n){if(!r[t]){if(!e[t]){var i="function"==typeof parcelRequire&&parcelRequire;if(!n&&i)return i(t,!0);if(o)return o(t,!0);if(u&&"string"==typeof t)return u(t);var c=new Error("Cannot find module '"+t+"'");throw c.code="MODULE_NOT_FOUND",c}p.resolve=function(r){return e[t][1][r]||r},p.cache={};var l=r[t]=new f.Module(t);e[t][0].call(l.exports,p,l,l.exports,this)}return r[t].exports;function p(e){return f(p.resolve(e))}}f.isParcelRequire=!0,f.Module=function(e){this.id=e,this.bundle=f,this.exports={}},f.modules=e,f.cache=r,f.parent=o,f.register=function(r,t){e[r]=[function(e,r){r.exports=t},{}]};for(var c=0;c<t.length;c++)try{f(t[c])}catch(e){i||(i=e)}if(t.length){var l=f(t[t.length-1]);"object"==typeof exports&&"undefined"!=typeof module?module.exports=l:"function"==typeof define&&define.amd?define(function(){return l}):n&&(this[n]=l)}if(parcelRequire=f,i)throw i;return f}({"K0KY":[function(require,module,exports) {
var define;
var t;!function(e,i,n,r){"use strict";var s,o=["","webkit","Moz","MS","ms","o"],a=i.createElement("div"),h="function",u=Math.round,c=Math.abs,l=Date.now;function p(t,e,i){return setTimeout(y(t,i),e)}function f(t,e,i){return!!Array.isArray(t)&&(v(t,i[e],i),!0)}function v(t,e,i){var n;if(t)if(t.forEach)t.forEach(e,i);else if(t.length!==r)for(n=0;n<t.length;)e.call(i,t[n],n,t),n++;else for(n in t)t.hasOwnProperty(n)&&e.call(i,t[n],n,t)}function d(t,i,n){var r="DEPRECATED METHOD: "+i+"\n"+n+" AT \n";return function(){var i=new Error("get-stack-trace"),n=i&&i.stack?i.stack.replace(/^[^\(]+?[\n$]/gm,"").replace(/^\s+at\s+/gm,"").replace(/^Object.<anonymous>\s*\(/gm,"{anonymous}()@"):"Unknown Stack Trace",s=e.console&&(e.console.warn||e.console.log);return s&&s.call(e.console,r,n),t.apply(this,arguments)}}s="function"!=typeof Object.assign?function(t){if(t===r||null===t)throw new TypeError("Cannot convert undefined or null to object");for(var e=Object(t),i=1;i<arguments.length;i++){var n=arguments[i];if(n!==r&&null!==n)for(var s in n)n.hasOwnProperty(s)&&(e[s]=n[s])}return e}:Object.assign;var m=d(function(t,e,i){for(var n=Object.keys(e),s=0;s<n.length;)(!i||i&&t[n[s]]===r)&&(t[n[s]]=e[n[s]]),s++;return t},"extend","Use `assign`."),g=d(function(t,e){return m(t,e,!0)},"merge","Use `assign`.");function T(t,e,i){var n,r=e.prototype;(n=t.prototype=Object.create(r)).constructor=t,n._super=r,i&&s(n,i)}function y(t,e){return function(){return t.apply(e,arguments)}}function E(t,e){return typeof t==h?t.apply(e&&e[0]||r,e):t}function I(t,e){return t===r?e:t}function A(t,e,i){v(b(e),function(e){t.addEventListener(e,i,!1)})}function _(t,e,i){v(b(e),function(e){t.removeEventListener(e,i,!1)})}function C(t,e){for(;t;){if(t==e)return!0;t=t.parentNode}return!1}function S(t,e){return t.indexOf(e)>-1}function b(t){return t.trim().split(/\s+/g)}function P(t,e,i){if(t.indexOf&&!i)return t.indexOf(e);for(var n=0;n<t.length;){if(i&&t[n][i]==e||!i&&t[n]===e)return n;n++}return-1}function D(t){return Array.prototype.slice.call(t,0)}function x(t,e,i){for(var n=[],r=[],s=0;s<t.length;){var o=e?t[s][e]:t[s];P(r,o)<0&&n.push(t[s]),r[s]=o,s++}return i&&(n=e?n.sort(function(t,i){return t[e]>i[e]}):n.sort()),n}function w(t,e){for(var i,n,s=e[0].toUpperCase()+e.slice(1),a=0;a<o.length;){if((n=(i=o[a])?i+s:e)in t)return n;a++}return r}var O=1;function R(t){var i=t.ownerDocument||t;return i.defaultView||i.parentWindow||e}var M="ontouchstart"in e,z=w(e,"PointerEvent")!==r,N=M&&/mobile|tablet|ip(ad|hone|od)|android/i.test(navigator.userAgent),X=25,Y=1,F=2,W=4,q=8,k=1,H=2,L=4,U=8,V=16,j=H|L,G=U|V,Z=j|G,B=["x","y"],$=["clientX","clientY"];function J(t,e){var i=this;this.manager=t,this.callback=e,this.element=t.element,this.target=t.options.inputTarget,this.domHandler=function(e){E(t.options.enable,[t])&&i.handler(e)},this.init()}function K(t,e,i){var n=i.pointers.length,s=i.changedPointers.length,o=e&Y&&n-s==0,a=e&(W|q)&&n-s==0;i.isFirst=!!o,i.isFinal=!!a,o&&(t.session={}),i.eventType=e,function(t,e){var i=t.session,n=e.pointers,s=n.length;i.firstInput||(i.firstInput=Q(e));s>1&&!i.firstMultiple?i.firstMultiple=Q(e):1===s&&(i.firstMultiple=!1);var o=i.firstInput,a=i.firstMultiple,h=a?a.center:o.center,u=e.center=tt(n);e.timeStamp=l(),e.deltaTime=e.timeStamp-o.timeStamp,e.angle=rt(h,u),e.distance=nt(h,u),function(t,e){var i=e.center,n=t.offsetDelta||{},r=t.prevDelta||{},s=t.prevInput||{};e.eventType!==Y&&s.eventType!==W||(r=t.prevDelta={x:s.deltaX||0,y:s.deltaY||0},n=t.offsetDelta={x:i.x,y:i.y});e.deltaX=r.x+(i.x-n.x),e.deltaY=r.y+(i.y-n.y)}(i,e),e.offsetDirection=it(e.deltaX,e.deltaY);var p=et(e.deltaTime,e.deltaX,e.deltaY);e.overallVelocityX=p.x,e.overallVelocityY=p.y,e.overallVelocity=c(p.x)>c(p.y)?p.x:p.y,e.scale=a?(f=a.pointers,v=n,nt(v[0],v[1],$)/nt(f[0],f[1],$)):1,e.rotation=a?function(t,e){return rt(e[1],e[0],$)+rt(t[1],t[0],$)}(a.pointers,n):0,e.maxPointers=i.prevInput?e.pointers.length>i.prevInput.maxPointers?e.pointers.length:i.prevInput.maxPointers:e.pointers.length,function(t,e){var i,n,s,o,a=t.lastInterval||e,h=e.timeStamp-a.timeStamp;if(e.eventType!=q&&(h>X||a.velocity===r)){var u=e.deltaX-a.deltaX,l=e.deltaY-a.deltaY,p=et(h,u,l);n=p.x,s=p.y,i=c(p.x)>c(p.y)?p.x:p.y,o=it(u,l),t.lastInterval=e}else i=a.velocity,n=a.velocityX,s=a.velocityY,o=a.direction;e.velocity=i,e.velocityX=n,e.velocityY=s,e.direction=o}(i,e);var f,v;var d=t.element;C(e.srcEvent.target,d)&&(d=e.srcEvent.target);e.target=d}(t,i),t.emit("hammer.input",i),t.recognize(i),t.session.prevInput=i}function Q(t){for(var e=[],i=0;i<t.pointers.length;)e[i]={clientX:u(t.pointers[i].clientX),clientY:u(t.pointers[i].clientY)},i++;return{timeStamp:l(),pointers:e,center:tt(e),deltaX:t.deltaX,deltaY:t.deltaY}}function tt(t){var e=t.length;if(1===e)return{x:u(t[0].clientX),y:u(t[0].clientY)};for(var i=0,n=0,r=0;r<e;)i+=t[r].clientX,n+=t[r].clientY,r++;return{x:u(i/e),y:u(n/e)}}function et(t,e,i){return{x:e/t||0,y:i/t||0}}function it(t,e){return t===e?k:c(t)>=c(e)?t<0?H:L:e<0?U:V}function nt(t,e,i){i||(i=B);var n=e[i[0]]-t[i[0]],r=e[i[1]]-t[i[1]];return Math.sqrt(n*n+r*r)}function rt(t,e,i){i||(i=B);var n=e[i[0]]-t[i[0]],r=e[i[1]]-t[i[1]];return 180*Math.atan2(r,n)/Math.PI}J.prototype={handler:function(){},init:function(){this.evEl&&A(this.element,this.evEl,this.domHandler),this.evTarget&&A(this.target,this.evTarget,this.domHandler),this.evWin&&A(R(this.element),this.evWin,this.domHandler)},destroy:function(){this.evEl&&_(this.element,this.evEl,this.domHandler),this.evTarget&&_(this.target,this.evTarget,this.domHandler),this.evWin&&_(R(this.element),this.evWin,this.domHandler)}};var st={mousedown:Y,mousemove:F,mouseup:W},ot="mousedown",at="mousemove mouseup";function ht(){this.evEl=ot,this.evWin=at,this.pressed=!1,J.apply(this,arguments)}T(ht,J,{handler:function(t){var e=st[t.type];e&Y&&0===t.button&&(this.pressed=!0),e&F&&1!==t.which&&(e=W),this.pressed&&(e&W&&(this.pressed=!1),this.callback(this.manager,e,{pointers:[t],changedPointers:[t],pointerType:"mouse",srcEvent:t}))}});var ut={pointerdown:Y,pointermove:F,pointerup:W,pointercancel:q,pointerout:q},ct={2:"touch",3:"pen",4:"mouse",5:"kinect"},lt="pointerdown",pt="pointermove pointerup pointercancel";function ft(){this.evEl=lt,this.evWin=pt,J.apply(this,arguments),this.store=this.manager.session.pointerEvents=[]}e.MSPointerEvent&&!e.PointerEvent&&(lt="MSPointerDown",pt="MSPointerMove MSPointerUp MSPointerCancel"),T(ft,J,{handler:function(t){var e=this.store,i=!1,n=t.type.toLowerCase().replace("ms",""),r=ut[n],s=ct[t.pointerType]||t.pointerType,o="touch"==s,a=P(e,t.pointerId,"pointerId");r&Y&&(0===t.button||o)?a<0&&(e.push(t),a=e.length-1):r&(W|q)&&(i=!0),a<0||(e[a]=t,this.callback(this.manager,r,{pointers:e,changedPointers:[t],pointerType:s,srcEvent:t}),i&&e.splice(a,1))}});var vt={touchstart:Y,touchmove:F,touchend:W,touchcancel:q},dt="touchstart",mt="touchstart touchmove touchend touchcancel";function gt(){this.evTarget=dt,this.evWin=mt,this.started=!1,J.apply(this,arguments)}T(gt,J,{handler:function(t){var e=vt[t.type];if(e===Y&&(this.started=!0),this.started){var i=function(t,e){var i=D(t.touches),n=D(t.changedTouches);e&(W|q)&&(i=x(i.concat(n),"identifier",!0));return[i,n]}.call(this,t,e);e&(W|q)&&i[0].length-i[1].length==0&&(this.started=!1),this.callback(this.manager,e,{pointers:i[0],changedPointers:i[1],pointerType:"touch",srcEvent:t})}}});var Tt={touchstart:Y,touchmove:F,touchend:W,touchcancel:q},yt="touchstart touchmove touchend touchcancel";function Et(){this.evTarget=yt,this.targetIds={},J.apply(this,arguments)}T(Et,J,{handler:function(t){var e=Tt[t.type],i=function(t,e){var i=D(t.touches),n=this.targetIds;if(e&(Y|F)&&1===i.length)return n[i[0].identifier]=!0,[i,i];var r,s,o=D(t.changedTouches),a=[],h=this.target;if(s=i.filter(function(t){return C(t.target,h)}),e===Y)for(r=0;r<s.length;)n[s[r].identifier]=!0,r++;r=0;for(;r<o.length;)n[o[r].identifier]&&a.push(o[r]),e&(W|q)&&delete n[o[r].identifier],r++;if(!a.length)return;return[x(s.concat(a),"identifier",!0),a]}.call(this,t,e);i&&this.callback(this.manager,e,{pointers:i[0],changedPointers:i[1],pointerType:"touch",srcEvent:t})}});var It=2500,At=25;function _t(){J.apply(this,arguments);var t=y(this.handler,this);this.touch=new Et(this.manager,t),this.mouse=new ht(this.manager,t),this.primaryTouch=null,this.lastTouches=[]}function Ct(t){var e=t.changedPointers[0];if(e.identifier===this.primaryTouch){var i={x:e.clientX,y:e.clientY};this.lastTouches.push(i);var n=this.lastTouches;setTimeout(function(){var t=n.indexOf(i);t>-1&&n.splice(t,1)},It)}}T(_t,J,{handler:function(t,e,i){var n="touch"==i.pointerType,r="mouse"==i.pointerType;if(!(r&&i.sourceCapabilities&&i.sourceCapabilities.firesTouchEvents)){if(n)(function(t,e){t&Y?(this.primaryTouch=e.changedPointers[0].identifier,Ct.call(this,e)):t&(W|q)&&Ct.call(this,e)}).call(this,e,i);else if(r&&function(t){for(var e=t.srcEvent.clientX,i=t.srcEvent.clientY,n=0;n<this.lastTouches.length;n++){var r=this.lastTouches[n],s=Math.abs(e-r.x),o=Math.abs(i-r.y);if(s<=At&&o<=At)return!0}return!1}.call(this,i))return;this.callback(t,e,i)}},destroy:function(){this.touch.destroy(),this.mouse.destroy()}});var St=w(a.style,"touchAction"),bt=St!==r,Pt="auto",Dt="manipulation",xt="none",wt="pan-x",Ot="pan-y",Rt=function(){if(!bt)return!1;var t={},i=e.CSS&&e.CSS.supports;return["auto","manipulation","pan-y","pan-x","pan-x pan-y","none"].forEach(function(n){t[n]=!i||e.CSS.supports("touch-action",n)}),t}();function Mt(t,e){this.manager=t,this.set(e)}Mt.prototype={set:function(t){"compute"==t&&(t=this.compute()),bt&&this.manager.element.style&&Rt[t]&&(this.manager.element.style[St]=t),this.actions=t.toLowerCase().trim()},update:function(){this.set(this.manager.options.touchAction)},compute:function(){var t=[];return v(this.manager.recognizers,function(e){E(e.options.enable,[e])&&(t=t.concat(e.getTouchAction()))}),function(t){if(S(t,xt))return xt;var e=S(t,wt),i=S(t,Ot);if(e&&i)return xt;if(e||i)return e?wt:Ot;if(S(t,Dt))return Dt;return Pt}(t.join(" "))},preventDefaults:function(t){var e=t.srcEvent,i=t.offsetDirection;if(this.manager.session.prevented)e.preventDefault();else{var n=this.actions,r=S(n,xt)&&!Rt[xt],s=S(n,Ot)&&!Rt[Ot],o=S(n,wt)&&!Rt[wt];if(r){var a=1===t.pointers.length,h=t.distance<2,u=t.deltaTime<250;if(a&&h&&u)return}if(!o||!s)return r||s&&i&j||o&&i&G?this.preventSrc(e):void 0}},preventSrc:function(t){this.manager.session.prevented=!0,t.preventDefault()}};var zt=1,Nt=2,Xt=4,Yt=8,Ft=Yt,Wt=16;function qt(t){this.options=s({},this.defaults,t||{}),this.id=O++,this.manager=null,this.options.enable=I(this.options.enable,!0),this.state=zt,this.simultaneous={},this.requireFail=[]}function kt(t){return t&Wt?"cancel":t&Yt?"end":t&Xt?"move":t&Nt?"start":""}function Ht(t){return t==V?"down":t==U?"up":t==H?"left":t==L?"right":""}function Lt(t,e){var i=e.manager;return i?i.get(t):t}function Ut(){qt.apply(this,arguments)}function Vt(){Ut.apply(this,arguments),this.pX=null,this.pY=null}function jt(){Ut.apply(this,arguments)}function Gt(){qt.apply(this,arguments),this._timer=null,this._input=null}function Zt(){Ut.apply(this,arguments)}function Bt(){Ut.apply(this,arguments)}function $t(){qt.apply(this,arguments),this.pTime=!1,this.pCenter=!1,this._timer=null,this._input=null,this.count=0}function Jt(t,e){return(e=e||{}).recognizers=I(e.recognizers,Jt.defaults.preset),new Kt(t,e)}qt.prototype={defaults:{},set:function(t){return s(this.options,t),this.manager&&this.manager.touchAction.update(),this},recognizeWith:function(t){if(f(t,"recognizeWith",this))return this;var e=this.simultaneous;return e[(t=Lt(t,this)).id]||(e[t.id]=t,t.recognizeWith(this)),this},dropRecognizeWith:function(t){return f(t,"dropRecognizeWith",this)?this:(t=Lt(t,this),delete this.simultaneous[t.id],this)},requireFailure:function(t){if(f(t,"requireFailure",this))return this;var e=this.requireFail;return-1===P(e,t=Lt(t,this))&&(e.push(t),t.requireFailure(this)),this},dropRequireFailure:function(t){if(f(t,"dropRequireFailure",this))return this;t=Lt(t,this);var e=P(this.requireFail,t);return e>-1&&this.requireFail.splice(e,1),this},hasRequireFailures:function(){return this.requireFail.length>0},canRecognizeWith:function(t){return!!this.simultaneous[t.id]},emit:function(t){var e=this,i=this.state;function n(i){e.manager.emit(i,t)}i<Yt&&n(e.options.event+kt(i)),n(e.options.event),t.additionalEvent&&n(t.additionalEvent),i>=Yt&&n(e.options.event+kt(i))},tryEmit:function(t){if(this.canEmit())return this.emit(t);this.state=32},canEmit:function(){for(var t=0;t<this.requireFail.length;){if(!(this.requireFail[t].state&(32|zt)))return!1;t++}return!0},recognize:function(t){var e=s({},t);if(!E(this.options.enable,[this,e]))return this.reset(),void(this.state=32);this.state&(Ft|Wt|32)&&(this.state=zt),this.state=this.process(e),this.state&(Nt|Xt|Yt|Wt)&&this.tryEmit(e)},process:function(t){},getTouchAction:function(){},reset:function(){}},T(Ut,qt,{defaults:{pointers:1},attrTest:function(t){var e=this.options.pointers;return 0===e||t.pointers.length===e},process:function(t){var e=this.state,i=t.eventType,n=e&(Nt|Xt),r=this.attrTest(t);return n&&(i&q||!r)?e|Wt:n||r?i&W?e|Yt:e&Nt?e|Xt:Nt:32}}),T(Vt,Ut,{defaults:{event:"pan",threshold:10,pointers:1,direction:Z},getTouchAction:function(){var t=this.options.direction,e=[];return t&j&&e.push(Ot),t&G&&e.push(wt),e},directionTest:function(t){var e=this.options,i=!0,n=t.distance,r=t.direction,s=t.deltaX,o=t.deltaY;return r&e.direction||(e.direction&j?(r=0===s?k:s<0?H:L,i=s!=this.pX,n=Math.abs(t.deltaX)):(r=0===o?k:o<0?U:V,i=o!=this.pY,n=Math.abs(t.deltaY))),t.direction=r,i&&n>e.threshold&&r&e.direction},attrTest:function(t){return Ut.prototype.attrTest.call(this,t)&&(this.state&Nt||!(this.state&Nt)&&this.directionTest(t))},emit:function(t){this.pX=t.deltaX,this.pY=t.deltaY;var e=Ht(t.direction);e&&(t.additionalEvent=this.options.event+e),this._super.emit.call(this,t)}}),T(jt,Ut,{defaults:{event:"pinch",threshold:0,pointers:2},getTouchAction:function(){return[xt]},attrTest:function(t){return this._super.attrTest.call(this,t)&&(Math.abs(t.scale-1)>this.options.threshold||this.state&Nt)},emit:function(t){if(1!==t.scale){var e=t.scale<1?"in":"out";t.additionalEvent=this.options.event+e}this._super.emit.call(this,t)}}),T(Gt,qt,{defaults:{event:"press",pointers:1,time:251,threshold:9},getTouchAction:function(){return[Pt]},process:function(t){var e=this.options,i=t.pointers.length===e.pointers,n=t.distance<e.threshold,r=t.deltaTime>e.time;if(this._input=t,!n||!i||t.eventType&(W|q)&&!r)this.reset();else if(t.eventType&Y)this.reset(),this._timer=p(function(){this.state=Ft,this.tryEmit()},e.time,this);else if(t.eventType&W)return Ft;return 32},reset:function(){clearTimeout(this._timer)},emit:function(t){this.state===Ft&&(t&&t.eventType&W?this.manager.emit(this.options.event+"up",t):(this._input.timeStamp=l(),this.manager.emit(this.options.event,this._input)))}}),T(Zt,Ut,{defaults:{event:"rotate",threshold:0,pointers:2},getTouchAction:function(){return[xt]},attrTest:function(t){return this._super.attrTest.call(this,t)&&(Math.abs(t.rotation)>this.options.threshold||this.state&Nt)}}),T(Bt,Ut,{defaults:{event:"swipe",threshold:10,velocity:.3,direction:j|G,pointers:1},getTouchAction:function(){return Vt.prototype.getTouchAction.call(this)},attrTest:function(t){var e,i=this.options.direction;return i&(j|G)?e=t.overallVelocity:i&j?e=t.overallVelocityX:i&G&&(e=t.overallVelocityY),this._super.attrTest.call(this,t)&&i&t.offsetDirection&&t.distance>this.options.threshold&&t.maxPointers==this.options.pointers&&c(e)>this.options.velocity&&t.eventType&W},emit:function(t){var e=Ht(t.offsetDirection);e&&this.manager.emit(this.options.event+e,t),this.manager.emit(this.options.event,t)}}),T($t,qt,{defaults:{event:"tap",pointers:1,taps:1,interval:300,time:250,threshold:9,posThreshold:10},getTouchAction:function(){return[Dt]},process:function(t){var e=this.options,i=t.pointers.length===e.pointers,n=t.distance<e.threshold,r=t.deltaTime<e.time;if(this.reset(),t.eventType&Y&&0===this.count)return this.failTimeout();if(n&&r&&i){if(t.eventType!=W)return this.failTimeout();var s=!this.pTime||t.timeStamp-this.pTime<e.interval,o=!this.pCenter||nt(this.pCenter,t.center)<e.posThreshold;if(this.pTime=t.timeStamp,this.pCenter=t.center,o&&s?this.count+=1:this.count=1,this._input=t,0===this.count%e.taps)return this.hasRequireFailures()?(this._timer=p(function(){this.state=Ft,this.tryEmit()},e.interval,this),Nt):Ft}return 32},failTimeout:function(){return this._timer=p(function(){this.state=32},this.options.interval,this),32},reset:function(){clearTimeout(this._timer)},emit:function(){this.state==Ft&&(this._input.tapCount=this.count,this.manager.emit(this.options.event,this._input))}}),Jt.VERSION="2.0.7",Jt.defaults={domEvents:!1,touchAction:"compute",enable:!0,inputTarget:null,inputClass:null,preset:[[Zt,{enable:!1}],[jt,{enable:!1},["rotate"]],[Bt,{direction:j}],[Vt,{direction:j},["swipe"]],[$t],[$t,{event:"doubletap",taps:2},["tap"]],[Gt]],cssProps:{userSelect:"none",touchSelect:"none",touchCallout:"none",contentZooming:"none",userDrag:"none",tapHighlightColor:"rgba(0,0,0,0)"}};function Kt(t,e){var i;this.options=s({},Jt.defaults,e||{}),this.options.inputTarget=this.options.inputTarget||t,this.handlers={},this.session={},this.recognizers=[],this.oldCssProps={},this.element=t,this.input=new((i=this).options.inputClass||(z?ft:N?Et:M?_t:ht))(i,K),this.touchAction=new Mt(this,this.options.touchAction),Qt(this,!0),v(this.options.recognizers,function(t){var e=this.add(new t[0](t[1]));t[2]&&e.recognizeWith(t[2]),t[3]&&e.requireFailure(t[3])},this)}function Qt(t,e){var i,n=t.element;n.style&&(v(t.options.cssProps,function(r,s){i=w(n.style,s),e?(t.oldCssProps[i]=n.style[i],n.style[i]=r):n.style[i]=t.oldCssProps[i]||""}),e||(t.oldCssProps={}))}Kt.prototype={set:function(t){return s(this.options,t),t.touchAction&&this.touchAction.update(),t.inputTarget&&(this.input.destroy(),this.input.target=t.inputTarget,this.input.init()),this},stop:function(t){this.session.stopped=t?2:1},recognize:function(t){var e=this.session;if(!e.stopped){var i;this.touchAction.preventDefaults(t);var n=this.recognizers,r=e.curRecognizer;(!r||r&&r.state&Ft)&&(r=e.curRecognizer=null);for(var s=0;s<n.length;)i=n[s],2===e.stopped||r&&i!=r&&!i.canRecognizeWith(r)?i.reset():i.recognize(t),!r&&i.state&(Nt|Xt|Yt)&&(r=e.curRecognizer=i),s++}},get:function(t){if(t instanceof qt)return t;for(var e=this.recognizers,i=0;i<e.length;i++)if(e[i].options.event==t)return e[i];return null},add:function(t){if(f(t,"add",this))return this;var e=this.get(t.options.event);return e&&this.remove(e),this.recognizers.push(t),t.manager=this,this.touchAction.update(),t},remove:function(t){if(f(t,"remove",this))return this;if(t=this.get(t)){var e=this.recognizers,i=P(e,t);-1!==i&&(e.splice(i,1),this.touchAction.update())}return this},on:function(t,e){if(t!==r&&e!==r){var i=this.handlers;return v(b(t),function(t){i[t]=i[t]||[],i[t].push(e)}),this}},off:function(t,e){if(t!==r){var i=this.handlers;return v(b(t),function(t){e?i[t]&&i[t].splice(P(i[t],e),1):delete i[t]}),this}},emit:function(t,e){this.options.domEvents&&function(t,e){var n=i.createEvent("Event");n.initEvent(t,!0,!0),n.gesture=e,e.target.dispatchEvent(n)}(t,e);var n=this.handlers[t]&&this.handlers[t].slice();if(n&&n.length){e.type=t,e.preventDefault=function(){e.srcEvent.preventDefault()};for(var r=0;r<n.length;)n[r](e),r++}},destroy:function(){this.element&&Qt(this,!1),this.handlers={},this.session={},this.input.destroy(),this.element=null}},s(Jt,{INPUT_START:Y,INPUT_MOVE:F,INPUT_END:W,INPUT_CANCEL:q,STATE_POSSIBLE:zt,STATE_BEGAN:Nt,STATE_CHANGED:Xt,STATE_ENDED:Yt,STATE_RECOGNIZED:Ft,STATE_CANCELLED:Wt,STATE_FAILED:32,DIRECTION_NONE:k,DIRECTION_LEFT:H,DIRECTION_RIGHT:L,DIRECTION_UP:U,DIRECTION_DOWN:V,DIRECTION_HORIZONTAL:j,DIRECTION_VERTICAL:G,DIRECTION_ALL:Z,Manager:Kt,Input:J,TouchAction:Mt,TouchInput:Et,MouseInput:ht,PointerEventInput:ft,TouchMouseInput:_t,SingleTouchInput:gt,Recognizer:qt,AttrRecognizer:Ut,Tap:$t,Pan:Vt,Swipe:Bt,Pinch:jt,Rotate:Zt,Press:Gt,on:A,off:_,each:v,merge:g,extend:m,assign:s,inherit:T,bindFn:y,prefixed:w}),(void 0!==e?e:"undefined"!=typeof self?self:{}).Hammer=Jt,"function"==typeof t&&t.amd?t(function(){return Jt}):"undefined"!=typeof module&&module.exports?module.exports=Jt:e.Hammer=Jt}(window,document);
},{}],"rqls":[function(require,module,exports) {
var global = arguments[3];
var t=arguments[3],e="undefined"!=typeof window?window:void 0!==t?t:"undefined"!=typeof self?self:{},n="Expected a function",i=NaN,r="[object Symbol]",o=/^\s+|\s+$/g,a=/^[-+]0x[0-9a-f]+$/i,s=/^0b[01]+$/i,u=/^0o[0-7]+$/i,l=parseInt,c="object"==typeof e&&e&&e.Object===Object&&e,f="object"==typeof self&&self&&self.Object===Object&&self,h=c||f||Function("return this")(),m=Object.prototype,d=m.toString,v=Math.max,p=Math.min,y=function(){return h.Date.now()};function g(t,e,i){var r,o,a,s,u,l,c=0,f=!1,h=!1,m=!0;if("function"!=typeof t)throw new TypeError(n);function d(e){var n=r,i=o;return r=o=void 0,c=e,s=t.apply(i,n)}function g(t){var n=t-l;return void 0===l||n>=e||n<0||h&&t-c>=a}function I(){var t=y();if(g(t))return C(t);u=setTimeout(I,function(t){var n=e-(t-l);return h?p(n,a-(t-c)):n}(t))}function C(t){return u=void 0,m&&r?d(t):(r=o=void 0,s)}function k(){var t=y(),n=g(t);if(r=arguments,o=this,l=t,n){if(void 0===u)return function(t){return c=t,u=setTimeout(I,e),f?d(t):s}(l);if(h)return u=setTimeout(I,e),d(l)}return void 0===u&&(u=setTimeout(I,e)),s}return e=w(e)||0,b(i)&&(f=!!i.leading,a=(h="maxWait"in i)?v(w(i.maxWait)||0,e):a,m="trailing"in i?!!i.trailing:m),k.cancel=function(){void 0!==u&&clearTimeout(u),c=0,r=l=o=u=void 0},k.flush=function(){return void 0===u?s:C(y())},k}function b(t){var e=typeof t;return!!t&&("object"==e||"function"==e)}function I(t){return!!t&&"object"==typeof t}function C(t){return"symbol"==typeof t||I(t)&&d.call(t)==r}function w(t){if("number"==typeof t)return t;if(C(t))return i;if(b(t)){var e="function"==typeof t.valueOf?t.valueOf():t;t=b(e)?e+"":e}if("string"!=typeof t)return 0===t?t:+t;t=t.replace(o,"");var n=s.test(t);return n||u.test(t)?l(t.slice(2),n?2:8):a.test(t)?i:+t}var k=g,N="function"==typeof Symbol&&"symbol"==typeof Symbol.iterator?function(t){return typeof t}:function(t){return t&&"function"==typeof Symbol&&t.constructor===Symbol&&t!==Symbol.prototype?"symbol":typeof t},E=function(t,e){if(!(t instanceof e))throw new TypeError("Cannot call a class as a function")},j=function(){function t(t,e){for(var n=0;n<e.length;n++){var i=e[n];i.enumerable=i.enumerable||!1,i.configurable=!0,"value"in i&&(i.writable=!0),Object.defineProperty(t,i.key,i)}}return function(e,n,i){return n&&t(e.prototype,n),i&&t(e,i),e}}(),A=function(t){if(Array.isArray(t)){for(var e=0,n=Array(t.length);e<t.length;e++)n[e]=t[e];return n}return Array.from(t)};function T(t){if(!t)return 0;var e=t.currentStyle||window.getComputedStyle(t),n=Number(e.width);return[Number.isNaN(n)?t.clientWidth:n,e.marginRight,e.marginLeft,e.paddingRight,e.paddingLeft,e.borderRightWidth,e.borderLeftWidth].reduce(function(t,e){var n=Number(e);return t+(Number.isNaN(n)?0:n)},0)}function W(t,e,n){var i=arguments.length>3&&void 0!==arguments[3]&&arguments[3];t.style.setProperty(e,n),i&&t.scrollTop}function S(t,e){var n=arguments.length>2&&void 0!==arguments[2]?arguments[2]:"div",i=document.createElement(n);return i.className=e,t.appendChild(i),i}function x(t){if("object"!==(void 0===t?"undefined":N(t))||!t.nodeType)throw new Error("Invalid element provided: "+(void 0===t?"undefined":N(t)));var e=[].concat(A(t.classList.values())).join(".");return t.tagName.toLowerCase()+(e&&e.length?".":"")+e}function D(t){if("object"===(void 0===t?"undefined":N(t))&&t.nodeType)return t;if("string"==typeof t)return document.querySelector(t);throw new Error("Invalid element or selector: "+t)}function q(t){var e=t.parentElement;return e?[].concat(A(e.querySelectorAll(x(t)+" > *"))):[]}var O={maxItems:800,duration:1e3,classNames:{container:"marqueeInfinite",slider:"marqueeInfiniteSlider",cell:"marqueeInfiniteCell"}},L=function(){function t(e,n){E(this,t);var i=D(e||".js-marquee-infinite");i&&(this.container=i,this.options=Object.assign({},O,n||{}),n&&n.classNames&&(this.options.classNames=Object.assign({},O.classNames,n.classNames)),this.setContainerClassName(),this.start())}return j(t,[{key:"setContainerClassName",value:function(){var t=this.options.classNames;this.container.classList.contains(t.container)||this.container.classList.add(t.container)}},{key:"start",value:function(){this.cacheContainerDirectChildren(),this.cacheInitialItems(),this.insertMarqueeElements(),this.cacheItems(),this.fillWithItems(),this.duplicateItems(),this.runAnimation(),this.events()}},{key:"refresh",value:function(){this.insertMarqueeElements(!0),this.cacheItems(),this.fillWithItems(),this.duplicateItems(),this.runAnimation()}},{key:"cacheContainerDirectChildren",value:function(){this.containerDirectChildren=q(this.container)}},{key:"cacheInitialItems",value:function(){if(this.initialItemElements=this.containerDirectChildren.map(function(t){return t.cloneNode(!0)}),!this.initialItemElements.length)throw new Error("Could not find any items");this.initialItemsWidth=this.containerDirectChildren.reduce(function(t,e){return t+T(e)},0)}},{key:"removeContainerDirectChildren",value:function(){this.containerDirectChildren.forEach(function(t){return t.parentElement.removeChild(t)})}},{key:"insertMarqueeElements",value:function(){var t=this;arguments.length>0&&void 0!==arguments[0]&&arguments[0]&&this.cacheContainerDirectChildren(),this.removeContainerDirectChildren();var e=this.options.classNames;this.slider=S(this.container,e.slider),this.container.appendChild(this.slider),this.initialItemElements.forEach(function(n){return S(t.slider,e.cell).appendChild(n)})}},{key:"cacheItems",value:function(){this.items=q(this.slider)}},{key:"addClone",value:function(t,e){e.appendChild(t.cloneNode(!0)),this.cacheItems()}},{key:"getAllItemsWidth",value:function(){return this.items.length*(this.initialItemsWidth/this.initialItemElements.length)}},{key:"getContainerWidth",value:function(){return T(this.container)}},{key:"getIsFilled",value:function(){return this.getAllItemsWidth()>=this.getContainerWidth()}},{key:"fillWithItems",value:function(){var t=this;if(!this.getIsFilled())for(var e=function(e){if(t.getIsFilled())return"break";var n=S(t.slider,t.options.classNames.cell);t.initialItemElements.forEach(function(e){return t.addClone(e,n)}),e>=t.options.maxItems-1&&!t.getIsFilled()&&console.warn("Reached the max items limit of "+t.options.maxItems+'. You can adjust this limit in your "options" object using the property "maxItems".')},n=0;n<this.options.maxItems;n+=this.initialItemElements.length){if("break"===e(n))break}}},{key:"duplicateItems",value:function(){var t=this;this.items.forEach(function(e){return t.addClone(e,t.slider)})}},{key:"getAnimationDuration",value:function(){return this.items.reduce(function(t,e){return t+T(e)},0)/(2*this.initialItemsWidth)*this.options.duration}},{key:"runAnimation",value:function(){var t=this;!function e(){var n=t.getAnimationDuration();W(t.slider,"transform","translateX(0%)"),W(t.slider,"transition","transform "+n+"ms linear"),W(t.slider,"transform","translateX(-100%)"),t.animationTimer=setTimeout(function(){t.stopAnimation(),requestAnimationFrame(function(){return e()})},n/2)}()}},{key:"stopAnimation",value:function(){W(this.slider,"transition","",!0),W(this.slider,"transform","translateX(0%)",!0)}},{key:"pauseAnimation",value:function(){this.animationTimer&&clearTimeout(this.animationTimer);var t=this.slider.currentStyle||window.getComputedStyle(this.slider);W(this.slider,"transform",t.transform),W(this.slider,"transition","",!0)}},{key:"events",value:function(){var t=this,e=this.handleResize.bind(this),n=k(e,1e3);window.addEventListener("resize",function(){t.pauseAnimation(),n()})}},{key:"handleResize",value:function(){this.refresh()}}]),t}();module.exports=L;
},{}],"Focm":[function(require,module,exports) {
"use strict";var e=t(require("hammerjs")),s=t(require("marquee-infinite"));function t(e){return e&&e.__esModule?e:{default:e}}new s.default(".m-container",{duration:6e4,classNames:{container:"marqueeInfinite",slider:"marqueeInfiniteSlider",cell:"marqueeInfiniteCell"}});var a=function(e){var s=document.getElementsByClassName("a")[0],t=document.getElementsByClassName("b")[0],a=document.getElementsByClassName("c")[0],m=(document.getElementsByClassName("d")[0],document.getElementsByClassName("e")[0],document.getElementsByClassName("f")[0]),l=document.getElementsByClassName("g")[0],n=document.getElementsByClassName("h")[0];s.classList.remove("a"),s.classList.add("b"),t.classList.remove("b"),t.classList.add("c"),a.classList.remove("c"),a.classList.add("f"),m.classList.remove("f"),m.classList.add("g"),l.classList.remove("g"),l.classList.add("h"),n.classList.remove("h"),n.classList.add("a")},m=function(e){var s=document.getElementsByClassName("a")[0],t=document.getElementsByClassName("b")[0],a=document.getElementsByClassName("c")[0],m=(document.getElementsByClassName("d")[0],document.getElementsByClassName("e")[0],document.getElementsByClassName("f")[0]),l=document.getElementsByClassName("g")[0],n=document.getElementsByClassName("h")[0];s.classList.remove("a"),s.classList.add("h"),t.classList.remove("b"),t.classList.add("a"),a.classList.remove("c"),a.classList.add("b"),m.classList.remove("f"),m.classList.add("c"),l.classList.remove("g"),l.classList.add("f"),n.classList.remove("h"),n.classList.add("g")},l=new e.default(document.getElementById("main"),{inputClass:e.default.TouchMouseInput});l.get("swipe").set({direction:e.default.DIRECTION_VERTICAL}),l.on("swipeup",a),l.on("swipedown",m);var n=function(e){38==e.keyCode?a():40==e.keyCode&&m()};document.addEventListener("keyup",n);
},{"hammerjs":"K0KY","marquee-infinite":"rqls"}]},{},["Focm"], null)
//# sourceMappingURL=/WorldsFirst.b317e33f.js.map