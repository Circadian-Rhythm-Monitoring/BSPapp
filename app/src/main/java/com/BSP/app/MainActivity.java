package com.BSP.app;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.util.Log;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;

import com.chaquo.python.PyObject;
import com.chaquo.python.Python;
import com.chaquo.python.android.AndroidPlatform;

public class MainActivity extends AppCompatActivity {

    PyObject output;
    PyObject output1;
    EditText nameInput;
    TextView helloworld;
    TextView outputDyn;
    TextView runtime;
    Button enter;

    public void outputMessage(View view){
        String name = nameInput.getText().toString();
        String msg = output.callAttr("helloWorld",name).toString();
        helloworld.setText(msg);
    }

    public void gainOptFilterDyn(View view){
        long startTime = System.nanoTime();
        Python py = Python.getInstance();

        Log.d("Python","Running Script");

        output1 = py.getModule("main");
        String test = output1.callAttr("main").toString();
        long endtime = System.nanoTime();

        outputDyn.setText("Output: " + test);
        Log.d("Python",test);
        long methodDuration = (endtime - startTime)/1_000_000_000;
        runtime.setText("RunTime: " + String.valueOf(methodDuration) + " seconds");
        Log.d("Python Runtime",String.valueOf(methodDuration));

    }


    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);


//        helloworld = findViewById(R.id.outputTextView);
//        nameInput = findViewById(R.id.editTextPersonName);
        enter = findViewById(R.id.buttonEnter);
        outputDyn = findViewById(R.id.textViewOutput);
        runtime = findViewById(R.id.textViewRuntime);

        if (! Python.isStarted()) {
            Python.start(new AndroidPlatform(this));
        }

//        Python py = Python.getInstance();
//        output = py.getModule("helloWorld");
//
//        Log.d("Python","Running Script");
//        output1 = py.getModule("main");
//        String test = output1.callAttr("main").toString();
//        Log.d("Python",test);


    }
}