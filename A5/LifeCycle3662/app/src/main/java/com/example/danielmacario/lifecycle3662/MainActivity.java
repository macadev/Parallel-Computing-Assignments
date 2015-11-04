package com.example.danielmacario.lifecycle3662;

import android.content.ComponentName;
import android.content.Context;
import android.content.Intent;
import android.content.ServiceConnection;
import android.os.Bundle;
import android.os.IBinder;
import android.support.v7.app.ActionBarActivity;
import android.util.Log;
import android.widget.Toast;

public class MainActivity extends ActionBarActivity {

    String msg = "Android : ";

    private Intent i;

    protected ServiceConnection mServerConn = new ServiceConnection() {
        @Override
        public void onServiceConnected(ComponentName name, IBinder binder) {
            Log.d("ServiceConn: ", "onServiceConnected");
        }

        @Override
        public void onServiceDisconnected(ComponentName name) {
            Log.d("ServiceConn: ", "onServiceDisconnected");
        }
    };

    /** Called when the activity is first created. */
    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        Log.d(msg, "The onCreate() event");

        displayToastMessage("The onCreate() event was called!");

        i = new Intent();
        i.setClassName("com.example.danielmacario.lifecycle3662",
                "com.example.danielmacario.lifecycle3662.MyService");
        bindService(i, mServerConn, Context.BIND_AUTO_CREATE);
        this.startService(i);

        setContentView(R.layout.activity_main);
    }

    /** Called when the activity is about to become visible. */
    @Override
    protected void onStart() {
        super.onStart();
        Log.d(msg, "The onStart() event");

        displayToastMessage("The onStart() event was called!");
    }

    /** Called when the activity has become visible. */
    @Override
    protected void onResume() {
        super.onResume();
        Log.d(msg, "The onResume() event");
        displayToastMessage("The onResume() event was called!");
    }

    /** Called when another activity is taking focus. */
    @Override
    protected void onPause() {
        super.onPause();
        Log.d(msg, "The onPause() event");
        displayToastMessage("The onPause() event was called!");
    }

    /** Called when the activity is no longer visible. */
    @Override
    protected void onStop() {
        super.onStop();
        Log.d(msg, "The onStop() event");
        displayToastMessage("The onStop() event was called!");
    }

    /** Called when the activity is restarted. */
    @Override
    protected void onRestart() {
        super.onRestart();
        Log.d(msg, "The onRestart() event");
        displayToastMessage("The onRestart() event was called!");
    }

    /** Called just before the activity is destroyed. */
    @Override
    public void onDestroy() {
        super.onDestroy();
        Log.d(msg, "The onDestroy() event");
        displayToastMessage("The onDestroy() event was called!");
        this.stopService(i);
        unbindService(mServerConn);
    }

    public void displayToastMessage(CharSequence notification) {
        Context context = getApplicationContext();
        int duration = Toast.LENGTH_SHORT;
        Toast toast = Toast.makeText(context, notification, duration);
        toast.show();
    }
}
