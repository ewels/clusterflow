//ZOMG! jQuery pong!! W00t, yeah.
//Based on the one by a guy named Ben White @ benwhite@columbus.rr.com
//jQuery'd by Ben Ogle. 

(function($){    
    $.fn.pong = function(ballImage, options) {
        
        var defaults = {
            targetSpeed: 30,    //ms
            ballAngle: 45,        //degrees
            ballSpeed: 8,         //pixels per update
            compSpeed: 5,         //speed of your opponent!!
            playerSpeed: 5,     //pixels per update
            difficulty: 5,
            width: 400,             //px
            height: 300,            //px
            paddleWidth: 10,    //px
            paddleHeight: 40, //px
            paddleBuffer: 1,    //px from the edge of the play area
            ballWidth: 14,        //px
            ballHeight: 14,     //px
            playTo: 10                //points
        }
        
        var opts = $.extend(defaults, options);
        
        ////functions!
        
        function PositionBall(bolSide, gameData, leftPaddle, rightPaddle, ball, score) {
            
            if (bolSide) {
                gameData.x = opts.width - opts.paddleWidth - opts.paddleBuffer - opts.ballWidth - 10;
            }
            else {
                gameData.x = opts.paddleWidth + opts.paddleBuffer + 10;
            }
            gameData.y = Math.round(Math.random() * (opts.height - ball.height()));
            
            ball.css('left', gameData.x);
            ball.css('top', gameData.y);
            
            if (bolSide != (0>Math.cos(opts.ballAngle*Math.PI/180)>0)) {
                opts.ballAngle += 180
            }
            
            ball.css('visibility', 'visible');
        }
        
        function UpdateScore(bolUser, gameData, leftPaddle, rightPaddle, ball, score, msg) {
            
            if (bolUser){
                gameData.playerScore++;
            }else{
                gameData.compScore++;
            }
            
            score.html('browser ' + gameData.compScore + ' | ' + 'you ' + gameData.playerScore);
            
            if (gameData.playerScore == opts.playTo || gameData.compScore == opts.playTo) {
                ball.css('visibility', 'hidden');
                gameData.gameOver = true;
                
                if(gameData.playerScore == opts.playTo)
                    score.append('; you win!');
                else
                    score.append('; you lose :(');
                    
            } else {
                PositionBall(bolUser, gameData, leftPaddle, rightPaddle, ball, score);
            }
        }
        
        ///Is run by the setTimeout function. Updates the gameData object. 
        function Update(gameData, leftPaddle, rightPaddle, ball, score, msg) {
            
            if (gameData.gameOver) {
                msg.html('click to start!');
                return;
            }
        
            msg.html('press ESC to stop');
            
            // Dynamically Adjust Game Speed
        
            var tmpDelay = new Date();
            var Diff = tmpDelay.valueOf() - gameData.delay.valueOf() - opts.target;
            gameData.speed += (Diff > 5)?-1:0;
            gameData.speed += (Diff < -5)?1:0;
            gameData.speed = Math.abs(gameData.speed);
            gameData.delay = tmpDelay;
        
            setTimeout(function(){Update(gameData, leftPaddle, rightPaddle, ball, score, msg)}, gameData.speed);
        
            //	MoveBall
        
            var d = opts.ballAngle * Math.PI / 180;
            gameData.y += Math.round(opts.ballSpeed*Math.sin(d));
            gameData.x += Math.round(opts.ballSpeed*Math.cos(d));
            var VB = 180-opts.ballAngle;
            var HB = 0-opts.ballAngle;
        
            //	Move Computer
        
            var LeftTop = parseInt(leftPaddle.css('top'));
            var LeftCenter = (opts.paddleHeight/2) + LeftTop
        
            if (Math.cos(d) > 0 || (gameData.x > opts.width/(2-(gameData.compAdj/(opts.difficulty*10))))) {
                var Center = (opts.height/2);
            } else {
                var BallTop = gameData.y;
                var Center = (opts.ballHeight/2) +BallTop;
            }
            var MoveDiff = Math.abs(Center - LeftCenter);
            if (MoveDiff > opts.compSpeed)
                MoveDiff = opts.compSpeed;
                
            if (Center > LeftCenter) 
                LeftTop += MoveDiff;
            else 
                LeftTop -= MoveDiff;

            if (LeftTop < 1)
                LeftTop = 1;
                
            if ((LeftTop+opts.paddleHeight+1) > opts.height) {
                LeftTop = opts.height - opts.paddleHeight - 1;
            }
            
            leftPaddle.css('top', LeftTop+'px');
        
            //	Move Player
        
            var RightTop = parseInt(rightPaddle.css('top'));
            if (gameData.up) 
                RightTop -= opts.playerSpeed;
            if (gameData.down) 
                RightTop += opts.playerSpeed;
    
            if (RightTop < 1)
                RightTop = 1;
            if ((RightTop+opts.paddleHeight+1) > opts.height) 
                RightTop=opts.height-opts.paddleHeight-1;

            rightPaddle.css('top', RightTop+'px');
        
            //	Check Top/Bottom/Left/Right
        
            if (gameData.y < 1) {
                gameData.y = 1;
                opts.ballAngle = HB;
            }
            
            if (gameData.y > opts.height-opts.ballHeight) {
                gameData.y = opts.height-opts.ballHeight;
                opts.ballAngle = HB;
            }
            
            if (gameData.x < 1) {
                gameData.x = 1;
                opts.ballAngle = VB;
                gameData.compAdj -= opts.difficulty;
                
                UpdateScore(true, gameData, leftPaddle, rightPaddle, ball, score, msg);
            }
            
            if (gameData.x > opts.width-opts.ballWidth) {
                gameData.x=opts.width-opts.ballWidth;
                opts.ballAngle=VB;
                UpdateScore(false, gameData, leftPaddle, rightPaddle, ball, score, msg);
            }
        
            //	Check Left Paddle
        
            var MaxLeft = opts.paddleWidth + opts.paddleBuffer;
            if (gameData.x < MaxLeft) {
                if (gameData.y < (opts.paddleHeight + LeftTop) && (gameData.y+opts.ballHeight) > LeftTop) {
                    gameData.x = MaxLeft;
                    opts.ballAngle = VB;
                    gameData.compAdj++;
                }
            }
        
            //	Check Right Paddle
        
            var MaxRight = opts.width - opts.ballWidth - opts.paddleWidth - opts.paddleBuffer;
            if (gameData.x > MaxRight) {
                if (gameData.y < (opts.paddleHeight + RightTop) && (gameData.y+opts.ballHeight) > RightTop) {
                    gameData.x = MaxRight;
                    opts.ballAngle = VB;
                }
            }
        
            ball.css('top', gameData.y);
            ball.css('left', gameData.x);
        
            if (gameData.compAdj < 0){
                gameData.compAdj = 0;
            }
        }
        
        function Start(gameData, leftPaddle, rightPaddle, ball, score, msg) {
            
            if (gameData.gameOver) {
                gameData.gameOver = false;
                gameData.playerScore = -1;
                gameData.compScore = -1;
                setTimeout(function(){Update(gameData, leftPaddle, rightPaddle, ball, score, msg)}, gameData.speed);
                UpdateScore(false, gameData, leftPaddle, rightPaddle, ball, score, msg);
                UpdateScore(true, gameData, leftPaddle, rightPaddle, ball, score, msg);
            }
        }
    
        return this.each(function() {
            
            var gameData = {
                up: false,                //Down key pressed?
                down: false,            //Down key pressed?
                x: 0,		//Ball X Pos
                y: 0,	                //Ball Y Pos
                compAdj: 0,	//Computer Adjust
                compScore: 0,	//Computer Score
                playerScore: 0,	//Player Score
                speed: 30,	//Actual Game Speed (Dynamic)
                gameOver: true,
                delay: new Date()
            }
            
            var $this = $(this);
            
            function keyDownEvent(event){
                switch (event.keyCode) {
                    case 38: //Up Arrow
                        gameData.up = true;
                        break;
                    case 40: //Down Arrow
                        gameData.down = true;
                        break;
                    case 27: //Escape
                        $this.children(".ball").css('visibility', 'hidden');
                        gameData.gameOver = true;
                        break;
                }
                return false;
            }
            
            function keyUpEvent(event){
                switch (event.keyCode) {
                    case 38: //Up Arrow
                        gameData.up = false;
                        break;
                    case 40: //Down Arrow
                        gameData.down = false;
                        break;
                }
                return false;
            }
            
            $this.css('background', '#000');
            $this.css('position', 'relative');
            
            $this.append('<textarea class="field" style="position:absolute;background:#000;border:0;top:-9999; left:-9999; width:0;height0;"></textarea>');
            $this.append('<div class="score" style="position:relative;color:#ffffff; font-family: sans-serif; text-align: center; font-weight: bold;">browser 0 | you 0</div>');
            $this.append('<div class="leftPaddle" style="position:absolute;background-color:#ffffff;"></div>');
            $this.append('<div class="rightPaddle" style="position:absolute;background-color:#ffffff;"></div>');
            $this.append('<img src="'+ballImage+'" class="ball" style="position:absolute;visibility:hidden;">');
            $this.append('<div class="msg" style="position:absolute; font-size: 8pt; color:#fff; bottom: 2px; right: 2px;"></div>');
            
            var leftPaddle = $this.children('.leftPaddle');
            var rightPaddle = $this.children('.rightPaddle');
            var ball = $this.children('.ball');
            var score = $this.children('.score');
            var msg = $this.children('.msg');
            var field = $this.children('.field');
            
            field.keydown( keyDownEvent );
            field.keyup( keyUpEvent );
            
            //field.css('width', 200);
            //field.css('height', 20);
            
            //initialize all
            $this.css('width', opts.width);
            $this.css('height', opts.height);

            leftPaddle.css('width', opts.paddleWidth);
            leftPaddle.css('height', opts.paddleHeight);
            leftPaddle.css('left', opts.paddleBuffer);
            leftPaddle.css('top', Math.round(1+(Math.random()*(opts.height-opts.paddleHeight-2))) );

            rightPaddle.css('width', opts.paddleWidth);
            rightPaddle.css('height', opts.paddleHeight);
            rightPaddle.css('left', opts.width - opts.paddleWidth - opts.paddleBuffer);
            rightPaddle.css('top', Math.round(1+(Math.random()*(opts.height-opts.paddleHeight-2))) );

            ball.css('width', opts.ballWidth);
            ball.css('height', opts.ballHeight);
            
            gameData.speed = opts.targetSpeed;
            Update(gameData, leftPaddle, rightPaddle, ball, score, msg);
            
            $this.click(function(){
                field.focus();
                Start(gameData, leftPaddle, rightPaddle, ball, score, msg);
            })
        });    
    };    
})(jQuery);    