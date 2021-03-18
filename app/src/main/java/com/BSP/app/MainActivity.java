package com.BSP.app;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;

import com.chaquo.python.PyObject;
import com.chaquo.python.Python;
import com.chaquo.python.android.AndroidPlatform;

public class MainActivity extends AppCompatActivity {

    PyObject output;
    EditText nameInput;
    TextView helloworld;
    Button enter;

    public void outputMessage(View view){
        String name = nameInput.getText().toString();
        String msg = output.callAttr("helloWorld",name).toString();
        helloworld.setText(msg);
    }


    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);


        helloworld = findViewById(R.id.outputTextView);
        nameInput = findViewById(R.id.editTextPersonName);
        enter = findViewById(R.id.buttonEnter);

        if (! Python.isStarted()) {
            Python.start(new AndroidPlatform(this));
        }

        Python py = Python.getInstance();
        output = py.getModule("helloWorld");









    }
}