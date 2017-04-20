from flask import Flask, render_template, request
import Gather_data

app = Flask(__name__)

@app.route("/")
def main():
    return render_template('index.html')
    
@app.route("/getvalues",methods = ['GET','POST'])
def getvalues():
	val = request.args.get('genename')
	val2 = Gather_data.get_names(val)
	#return 'hello'
	return val2
    

if __name__ == "__main__":
    app.run()
