// this is the only part that gets used.

var bufsize = 2048;
var ping = [];
var ping_length = 3;
var global_elapsed;
var buftime = bufsize / 44100;

if (Pizzicato.Sound) { 
    var ticker = new Pizzicato.Sound({
        source: 'script',
        options: {
            attack: 0,
            release: 0.1,
            bufferSize: bufsize,
            volume: 0.2,
            audioFunction: function(e) {
                output = e.outputBuffer.getChannelData(0);
                var l = e.outputBuffer.length;
                // flat to start...
                for (var i = 0; i < l; i++) {
                  output[i] = -1;
                }
                ping.sort();
                var el = global_elapsed;
                var bufend = el + buftime;
                var next_ping = ping.shift();
                while (next_ping && next_ping < bufend) {
                  var st = Math.floor((el - next_ping)*44.100);
                  for (var i=0; i<ping_length && st >=0 && (st+i)<bufsize; i++) {
                    output[st+i] = 1;
                  }
                  next_ping = ping.shift();
                }
                // put it back if there's one left over.
                if (next_ping) { ping.splice(0, 0, next_ping) }
            }
        }
    });
} else {
    var ticker = null;
}
  
