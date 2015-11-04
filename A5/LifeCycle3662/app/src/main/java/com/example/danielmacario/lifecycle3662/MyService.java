package com.example.danielmacario.lifecycle3662;

import android.app.Service;
import android.content.Context;
import android.content.Intent;
import android.database.ContentObserver;
import android.os.Handler;
import android.os.IBinder;
import android.util.Log;
import android.widget.Toast;

public class MyService extends Service {

    private String tag = "ServiceNotification: ";
    private ContentObserver mSettingsContentObserver;

    public MyService() {
    }

    Handler volumeChangeHandler;

    @Override
    public void onCreate() {
        super.onCreate();
        Log.i(tag, "Service created.");
        displayToastMessage("Volume Tracker Service created!");

        volumeChangeHandler = new Handler();

        mSettingsContentObserver = new SettingsContentObserver(this, volumeChangeHandler);
        getApplicationContext().getContentResolver().registerContentObserver(android.provider.Settings.System.CONTENT_URI, true, mSettingsContentObserver);
    }

    @Override
    public void onDestroy() {
        super.onDestroy();
        Log.i(tag, "Service destroyed");
        displayToastMessage("Volume Tracker Service destroyed!");

        getApplicationContext().getContentResolver().unregisterContentObserver(mSettingsContentObserver);

    }

    public void displayToastMessage(CharSequence notification) {
        Context context = getApplicationContext();
        int duration = Toast.LENGTH_SHORT;
        Toast toast = Toast.makeText(context, notification, duration);
        toast.show();
    }

    @Override
    public IBinder onBind(Intent intent) {
        // TODO: Return the communication channel to the service.
        //throw new UnsupportedOperationException("Not yet implemented");
        return null;
    }
}
